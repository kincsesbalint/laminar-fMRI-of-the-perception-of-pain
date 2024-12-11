function [columnwisestat,columndistribution,layerwisestat]=BK_layer_sampling_pain_study_pipeline(subid)
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




% These subject are ready:
% 7349,                    7361,  7375,  7376(!!!),  7383,  7402,  7403,  7404,  7405
%
% List of all subjects:
% 7349,  7356(!),  7361,  7375,  7376(!!!),  7383,  7402,  7403,  7404,  7405
% 7408,  7415,  7426,  7434,  7443,  7445,  7449,  7453,  7455,  7457,  7469,  7485  
% 7414,  7425,  7433,  7435,  
% 7444,  7448,  7452,  7454,  7456,  7468,  7482,  7484
% 
%7356the orientation is incorrect, therefore the matlab sapce htingy is not
%good....
%
% 
% Example call of the function: diary subj8subj_layerif; for subj={7408,  7415,  7426,  7434,  7443,  7445,  7449,  7453}; BK_layer_sampling_pain_study_pipeline(subj); end; diary off

fspath = 'E:\pain_layers\main_project\derivatives\freesurfer\';

tic
% for subid = subids
%     subid
    if subid{1} <= 7405
        subpath = 'E:/pain_layers/main_project/derivatives/pipeline/';
    else
        subpath = 'D:/main_project/derivatives/pipeline/';
    end
    if subid{1}==7356
        T1path = [subpath num2str(subid{1}) '' ...
              '\ses-01\anat\presurf_MPRAGEise\presurf_UNI\UNI_MPRAGEised_biascorrected.nii'];
        
    else
        T1path = [subpath num2str(subid{1}) '' ...
              '\ses-01\anat\presurf_MPRAGEise\presurf_UNI\UNI_MoCo_MPRAGEised_biascorrected.nii'];
    end
        BK_checkfile(T1path)
%     BK_ROI_creation_GROUPLVLactivationconj(subid{:},subpath,fspath,T1path,0,'rwls'); %no visualization
    [columnwisestat,columndistribution,layerwisestat]=BK_firstlvlanalysis(subid{:},subpath,fspath,T1path,0,'rwls'); %no visualization

%     BK_ROI_creation_GROUPLVLactivationconj(subid{:},subpath,fspath,T1path,1); %with visualization
    fprintf('Time for one subject to process:')
    toc
end

