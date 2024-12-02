% This is a modified script of VPF's same function. 
% I removed many part of it. Thos parts which has saved the output in a
% permanent location (in the external HDs), were removed. Additionally,
% output warning messages when these files are somehow not available. 
% 
% The function calls
% 
% Balint Kincses
% 08.2024
% 

clear;clc;
%32+1=33 subjects are run. except subj 7376
subids ={};
%

%only the ROI creation was run for these subjects:

%the whole pipeline was run for these subjects:
% subids = {'7408','7414','7415'};
% subids ={'7425','7426','7433','7434','7435','7443','7444','7445','7448','7449', '7452','7453'};
% subids ={'7454','7455','7456','7457','7468','7469'};
% subids ={'7482','7484','7349','7361','7375'};
% subids ={'7383','7402','7403','7404','7405','7356','7485'};

fspath = 'E:\pain_layers\main_project\derivatives\freesurfer\';
tic
for subid = subids
    subid
    if str2double(subid{:}) <= 7405
        subpath = 'E:/pain_layers/main_project/derivatives/pipeline/';
    else
        subpath = 'D:/main_project/derivatives/pipeline/';
    end
    if str2double(subid{:})==7356
        T1path = [subpath subid{:} '' ...
              '\ses-01\anat\presurf_MPRAGEise\presurf_UNI\UNI_MPRAGEised_biascorrected.nii'];
        
    else
        T1path = [subpath subid{:} '' ...
              '\ses-01\anat\presurf_MPRAGEise\presurf_UNI\UNI_MoCo_MPRAGEised_biascorrected.nii'];
    end
        BK_checkfile(T1path)
    
    % Elapsed time is 65.837066 seconds for ROI creation:
    %Elapsed time is 73.456412 seconds.

%     BK_ROI_creation(subid{:},subpath,fspath);

    %37.9366min for one subject
%     tic
    layers = BK_layer_sampling_pain_study(subid{:},subpath,fspath,T1path);
%     display('Time for layer sampling:')
%     toc
% %     
    save([subpath '/' subid{:} '/ses-02/func/layers/smoothed_layers.mat'],'layers');
end
toc
