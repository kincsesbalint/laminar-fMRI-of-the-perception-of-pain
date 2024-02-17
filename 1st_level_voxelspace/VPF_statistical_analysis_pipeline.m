%this script performes robust weighted linear least squares analysis on 
%multiple subjects and estimates the correspoding SPM.mat.

%The output will be a folder named rwls_stats for each experiment (WM_localizer,
%pain_calib, post_calib, layers) and subject. In this folder, the SPM.mat and
%beta maps are saved.
clear;clc;

pipepath = '/media/pfaffenrot/My Passport2/main_project/derivatives/pipeline/';

subjects = 7469;
experiments = [cellstr('WM_localizer'),cellstr('pain_calib'),cellstr('post_calib')];
% cellstr('WM_localizer'),cellstr('pain_calib'),cellstr('post_calib')   
for subject = subjects
    for experiment = experiments
        VPF_execute_SPM_batch(pipepath,subject,experiment{:});
    
        %compress the used nifies to save space
        experiment_path = [pipepath num2str(subject) '/ses-02/func/' experiment{:}];
        Nruns = length(dir([experiment_path '/run*']));

        for run = 1:Nruns
            dat_path = [experiment_path '/run' num2str(run) '/func/'];
            files = dir([dat_path '/*Warped-to-Anat.nii']);

            for ii = 1:length(files)
                file = [files(ii).folder '/' files(ii).name];
                system(['pigz ' '"' file '"']);
            end

        end
    end
end


