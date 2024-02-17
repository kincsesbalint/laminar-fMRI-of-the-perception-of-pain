clear;clc;

subpath = '/media/pfaffenrot/My Passport2/main_project/derivatives/pipeline';
fspath = '/media/pfaffenrot/My Passport1/pain_layers/main_project/derivatives/freesurfer';

subids = {'7453','7454','7455','7456','7457'};




for subid = subids
T1path = ['/media/pfaffenrot/My Passport2/main_project/derivatives/pipeline/' subid{:} '' ...
          '/ses-01/anat/presurf_MPRAGEise/presurf_UNI/UNI_MoCo_MPRAGEised_biascorrected.nii'];


VPF_ROI_creation(subid{:},subpath,fspath);
layers = VPF_layer_sampling_pain_study(subid{:},subpath,fspath,T1path);

save([subpath '/' subid{:} '/ses-02/func/layers/layers.mat'],'layers');
end
