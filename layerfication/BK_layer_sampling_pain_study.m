function layers = BK_layer_sampling_pain_study(subid,subpath,fspath,T1path,N_layers)
% This is a adaptation of VPF's similarly called function.
%Function to sample layers within predefined regions (see VPF_ROI_creation.m). 
%It loads freesurfer surfaces, transforms ROI volumetric .nii to surfaces 
% and samples layers of images saved in folders with 'run' in their name
%INPUT:
%subid [str]                : subject id
%subpath [str]              : path to the derivatives folder
%fspath [str]               : path to the freesurfer folder
%T1path [str]               : full path and filename of the T1 image
%N_layers [int] optional    :number of layers to be sampled (default: 20)
%OUTPUT:
%layer [N_layers N_ROI N_vol N_run] : sampled layer array
%
% 
% 
% Balint Kincses
% 08.2024
if ~ischar(subid)
    subid = char(string(subid));
end
if nargin < 5
    N_layers = 20;
end
tic
%load boundaries coming from freesurfer
layer_boundaries = VPF_load_layer_boundaries(subid,fspath);
fprintf('time for loading:')
toc
%load ROIs coming from "VPF_ROI_creation.m"
% time for converting:Elapsed time is 92.710407 seconds.
tic
ROIs = VPF_convert_load_ROIs(subid,subpath,fspath,size(layer_boundaries,1));
fprintf('time for converting:')
toc
%transform boundaries to matlab space
tic
fprintf('time for transforming:')
layer_boundaries = VPF_transform_layers_to_matlab_space(layer_boundaries,T1path);
toc
%sample layers
%todo check the code from here:
fprintf(sprintf('Starting layer sampling...\n'));
layers = VPF_sample_layers(subid,subpath,layer_boundaries,ROIs,N_layers);
end

function layer_boundaries = VPF_load_layer_boundaries(subid,fspath)
%loads the layer boundaries (white/pial) generated by Freesurfer
%INPUT:
%subid [str]                : subject id
%fspath [str]               : path to the freesurfer folder
%OUTPUT:
% layer_boundaries [vertex,[x,y,z],[white,pial]] : 3D matrix of the surfaces
%this part loads in the vertex_coordinates of the surface file
%(https://github.com/fieldtrip/fieldtrip/blob/master/external/freesurfer/read_surf.m)
%and create a huge matrix which each row represent the cooridate of a
%vertex (starting from left hemisphere and continou with right hemishpere)
%and teh 3rd dimension contains white matter and pial information.
fn = [fspath '\' subid '\surf\lh.white'];
layer_boundaries = read_surf(fn);
fn = [fspath '\' subid '\surf\rh.white'];
tmp = read_surf(fn);
layer_boundaries  = cat(1,layer_boundaries ,tmp);


fn = [fspath '\' subid '\surf\lh.pial'];
tmp = read_surf(fn);
fn = [fspath '\' subid '\surf\rh.pial'];
tmp = cat(1,tmp,read_surf(fn));
layer_boundaries = cat(3,layer_boundaries,tmp);

end

function ROIs = VPF_convert_load_ROIs(subid,subpath,fspath,N_vertex)
%converts the ROIs produced by "VPF_ROI_creation.m" to surfaces and loads them
%INPUT:
%subid [str]                : subject id
%subpath [str]              : path to the derivatives folder
%fspath [str]               : path to the freesurfer folder


%OUTPUT:
%ROIs [N_vertex,N_ROI] : 2D logical matrix containing the ROIs in surface
%                        space. ROIs are sorted the following way: First
%                        come all the WM localizer, then the pain localizer.
%                        Within each type, the ROIS are sorted by ROI number
%                        in ascending order. Each ROI is split into voxels
%                        with a negative response and a positive response.
%                        The final structure is WM_ROI1_neg,WM_ROI1_pos,
%                        WM_ROI2_neg,....,pain_ROI1_neg,pain_ROI1_pos,...

% fs_base = ['export FREESURFER_HOME=/usr/local/freesurfer/6.0.0; ' ...
%     'export SUBJECTS_DIR=' fspath '; '...
%     'source $FREESURFER_HOME/SetUpFreeSurfer.sh; '];


%working memory localizer
WM_localizer = dir(fullfile([subpath '\' subid '\ses-02\func\layers\rwls_stats_compcor_smoothed'], '**', 'WM_localizer_ROI*.nii.gz'));

%pain localizer
pain_localizer = dir(fullfile([subpath '\' subid '\ses-02\func\layers\rwls_stats_compcor_smoothed'], '**', 'pain_localizer_ROI*.nii.gz'));

N_total_ROIs = length(WM_localizer)+length(pain_localizer);
localizer = cat(1,WM_localizer,pain_localizer);

%first convert freesurfer's xh.white to xh.white.gii
%theoretically this step is needed so wb_command then capable of using the
%surface data (othe solruion would be to use mri_vol2surf??) but most
%probably that result in different result
%(https://www.mail-archive.com/freesurfer@nmr.mgh.harvard.edu/msg66816.html)

%check if these files exist
% cmd = ['mris_convert ' '"' fspath '/' subid '/surf/lh.white" lh.white.gii'];
% system([fs_base cmd]);
BK_checkfile([fspath '\' subid '\surf\lh.white.gii'])
% cmd = ['mris_convert ' '"' fspath '/' subid '/surf/rh.white" rh.white.gii'];
% system([fs_base cmd]);
BK_checkfile([fspath '\' subid '\surf\rh.white.gii'])

ROIs = zeros([N_vertex,N_total_ROIs],'logical');
for ROI = 1:N_total_ROIs
    data = [localizer(ROI).folder '\' localizer(ROI).name];
    wls_data=strrep(data, '\', '/');
    if data(1)=='E'
        wls_data=strrep(wls_data,'E:','/mnt/e');
    elseif data(1)=='D'
        wls_data=strrep(wls_data,'D:','/mnt/d');
    end
    wls_fspath=strrep(fspath, '\', '/');
    wls_fspath=strrep(wls_fspath,'E:','/mnt/e');
    
    try %try to load to ROI surfaces
        ROIs(:,ROI) = VPF_average_mask_and_thres(data(1:end-4));
    catch %apparently, not created yet so create them first
        hemis = {'lh','rh'};
        for hemi = 1:2
            cmd = ['wsl wb_command -volume-to-surface-mapping ' '"' wls_data '"'...
                ' "' wls_fspath '/' subid '/surf/' hemis{hemi} '.white.gii' '"'...
                ' "' wls_data(1:end-4) '_' hemis{hemi} '.shape.gii' '"'...
                ' -cubic'];
            system(cmd);
        end
        ROIs(:,ROI) = VPF_average_mask_and_thres(data(1:end-4));
    end
end
end

function cdata = VPF_average_mask_and_thres(in)
%loads the ROI surfaces. As they are transformed from voxel to surface space, 
%some interpolation errors occur which are fixed here. Finally, the left 
%and right hemispheres are concatenated
%INPUT:
%in [str]                : full path and name of ROI WITHOUT extension (i.e.
%                          without .nii)      


%OUTPUT:
%cdata [N_vertex]        : logical array of the ROI in surface space

ROI = gifti([in '_lh.shape.gii']);
cdata_lh = ROI.cdata;
cdata_lh(abs(cdata_lh)>=1e-2) = 1;
cdata_lh(abs(cdata_lh)<1e-2) = 0;

ROI = gifti([in '_rh.shape.gii']);
cdata_rh = ROI.cdata;
cdata_rh(abs(cdata_rh)>=1e-2) = 1;
cdata_rh(abs(cdata_rh)<1e-2) = 0;

cdata = logical(cat(1,cdata_lh,cdata_rh));
end

function layer_boundaries = VPF_transform_layers_to_matlab_space(layer_boundaries,T1path)
%Transforms the layers from freesurfer space to matlab space. After that, 
%signal can be sampled using spm_sample_vol
%INPUT:
%layer_boundaries [vertex,[x,y,z],[white,pial]] : 3D matrix of the surfaces
%T1path [str]                                   : full path and filename of 
%                                                 the T1 image

%OUTPUT:
%layer_boundaries        : 3D matrix of surfaces in matlab space

sz = size(layer_boundaries);

%bring surfaces into matlab space
T1_mat = spm_get_space(T1path);
boundaries_2_cat = ones([sz(1) 1 sz(end)]);
layer_boundaries = cat(2,layer_boundaries,boundaries_2_cat);

for ii = 1:sz(end)
    tmp = T1_mat\squeeze(layer_boundaries(:,:,ii)).';
    layer_boundaries(:,:,ii) = tmp.';
end

layer_boundaries= layer_boundaries(:,1:3,:);

% img_T1 = spm_read_vols(spm_vol(T1path));
% for slice = 205
%     idx = find(layer_boundaries(:,3,1)>slice & layer_boundaries(:,3,1)<slice+1);
%     idx3 = find(layer_boundaries(:,3,2)>slice & layer_boundaries(:,3,2)<slice+1);
% 
%     figure;imagesc(squeeze(img_T1(:,:,slice)), [0 1200]), colormap(gray(256)), title(['Slice ' num2str(slice)])
%     hold on;
% 
%     plot(layer_boundaries(idx,2,1),layer_boundaries(idx,1,1),'r.');
%     plot(layer_boundaries(idx3,2,2),layer_boundaries(idx3,1,2),'y.');
% 
%     legend('white','pial');
% end

end

function layers_out = VPF_sample_layers(subid,subpath,layer_boundaries,ROIs,N_layers)
%Core sampling function. Searches for folders with 'run' in their name and
%samples the containing *Warped-to-Anat.nii.gz images within the predefined 
%ROIs.

%INPUT:
%subid [str]                : subject id
%subpath [str]              : path to the derivatives folder
%layer_boundaries           : 3D matrix of surfaces in matlab space
%ROIs [N_vertex,N_ROI] :    : 2D logical matrix containing the ROIs in
%                             surface space
%N_layers [int]             : number of layers (default: 20)

%OUTPUT:
%layer_out [N_layers N_ROI N_vol N_run]
%                           : sampled layer array
currpath = pwd;
SAMPLING_ORDER = 3;
runs = dir([subpath '\' subid '\ses-02\func\layers\run*']);
N_run = numel(runs);
N_ROI = size(ROIs,2);
for run = 1:N_run
    fprintf(sprintf(['run ' num2str(run) '...']));

    cd([runs(run).folder '\' runs(run).name '\func']);
    sampled_img_list = dir([runs(run).folder '\' runs(run).name '\func\*Warped-to-Anat.nii.gz']);
    N_vols = length(sampled_img_list);

    if run==1
        layers_out = zeros([N_layers,N_ROI,N_vols,N_run]);
    end
    for vol = 1:N_vols
%         sampled_img = load_nifti([sampled_img_list(vol).name]).vol;
        sampled_img = niftiread([sampled_img_list(vol).name]);
        parfor ROI = 1:N_ROI
            ind = find(ROIs(:,ROI));
            tmp = zeros(size(ind,1),N_layers);
            boundaries = layer_boundaries;
%we go here from vertex to vertex(==column to column), and sample the funcitonal images. Here we should include an additional function, to select only columns active mostly in the deep layer. However, this need to be happen on the unsmoothed main effects (smoothed actually is available) conjuction map. That is calculate indiviudal conjuction maps, and apply the mask on those and estimate the layer profile. However, we would need to have a decent t value at the deeper part of the cortex. (my gut feeleing that interpolate statistical map is not that correct.-think about it a bit more)
            for k = 1:size(ind,1) 
                X  = linspace(boundaries(ind(k),1,1),boundaries(ind(k),1,2),N_layers);
                Y  = linspace(boundaries(ind(k),2,1),boundaries(ind(k),2,2),N_layers);
                Z  = linspace(boundaries(ind(k),3,1),boundaries(ind(k),3,2),N_layers);
                
                tmp(k,:) = spm_sample_vol(sampled_img,X,Y,Z,SAMPLING_ORDER);
            end
            layers_out(:,ROI,vol,run) = squeeze(mean(tmp,1,'omitnan'));
        end
    end
    fprintf(sprintf(' done \n'));
end
cd(currpath)
end

