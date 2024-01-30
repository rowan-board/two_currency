
%% Script to make time series file over region of interest (txt file)

% Note in order for the fsl commands to work, MATLAB needs to be launched
% from terminal e.g. /Applications/MATLAB_R2021a.app/bin/matlab &

%% set FSL environment
setenv('FSLDIR','/usr/local/fsl/');  % location of FSL folder
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % type of output 

%% Define files

subjList = [ 03 04 ]; 
ROI_NAME = 'vmPFC'; % name of mask to run time series on 

dataFolder = '/Users/rowanboard/Desktop/analyse_timeCourses/'; 
mask_directory = '/Users/rowanboard/Desktop/analyse_timeCourses/';

maskFolder = strcat([mask_directory ROI_NAME]);

standard2highres_warps_folder = [mask_directory,'standard2highres_warps'];
if ~isfolder(standard2highres_warps_folder)
    mkdir(standard2highres_warps_folder);
end

individual_masks_folder = [maskFolder, '/individual_masks'];
if ~isfolder(individual_masks_folder)
    mkdir(individual_masks_folder);
end

timeCourse_files_folder = [maskFolder, '/timeCourse_files'];
if ~isfolder(timeCourse_files_folder)
    mkdir(timeCourse_files_folder);
end 

standard_space_mask = [maskFolder, '/',ROI_NAME,'_sphere.nii.gz'];
%% Loop over subjects and runs

for s = 1:numel(subjList) 
    
    for run = 1:2
    
        ID = subjList(s); 

        % Individual registration data path
        subPath = [dataFolder 'sub_0',num2str(ID), '/run_', num2str(run)];

        % Denoised functional file 
        func_denoised = [subPath '/filtered_func_data.nii.gz'];

        % Registration files 
        highres2standard_warp = [subPath '/reg/highres2standard_warp.nii.gz'];
        highres = [subPath '/reg/highres.nii.gz'];
        example_func = [subPath '/reg/example_func.nii.gz'];
        highres2example_func = [subPath '/reg/highres2example_func.mat'];

        % Name of warp file to go from standard space to subject space
        standard2highres_warp = [standard2highres_warps_folder, '/sub_0' num2str(ID) '_run' num2str(run) '_standard2highres.nii.gz']; 

        % Create standard —> subject space warp file if it doesn't already exist
        if ~isfile(standard2highres_warp)
            cmdStr1 = ['/usr/local/fsl/bin/invwarp -w ' highres2standard_warp ' -o ' standard2highres_warp ' -r ' highres];
            system(cmdStr1)
        end

        individual_mask_file = [individual_masks_folder '/sub_0' num2str(ID) '_run' num2str(run)]; 

        % Use the new registration file to create the individual subject mask 
        cmdStr2 = ['/usr/local/fsl/bin/applywarp -i ' standard_space_mask ' -r ' example_func ' -o ' individual_mask_file ' -w ' standard2highres_warp ' --postmat=' highres2example_func];
        system(cmdStr2) 

        individual_mask_file_thresholded = [individual_mask_file '_bin'];

        % Threshold the mask - note to ask what would be a good value to
        % threshold individual masks to? 
%         cmdStr3 = ['fslmaths ' individual_mask_file ' -thr 0.5 -bin ' individual_mask_file_thresholded];
%         system(cmdStr3)
        
        time_course_file = [timeCourse_files_folder '/sub_0' num2str(ID) '_run' num2str(run)];
        
        % Use the mask to create a textfile of the extracted (and averaged) timecourse from the ROI
        cmdStr4 = ['/usr/local/fsl/bin/fslmeants -i ' func_denoised ' -m ' individual_mask_file ' -o ' time_course_file '.txt'];
        %cmdStr4 = ['fslmeants -i ' func_denoised ' -m ' individual_mask_file_thresholded ' -o ' individual_mask_file '.txt'];
        system(cmdStr4)
    end
    
    
end


    