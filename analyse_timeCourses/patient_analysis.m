%% patient vs control analyses 
%% Define paths 

homepath = '/Users/Rowan/Desktop/oxford';
addpath(genpath(homepath)); 

neural_dat = [homepath '/analyse_timeCourses/'];

% voxel wise (1) or mean ROI (0)?
voxel = 0;

participants = {'02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23'}; 

if voxel==0
    ROIs = {'vmPFC','dmPFC','VS','mFP','lFP'};
    ROIs = {'bilat_vmpfc_roi', 'bilat_vs_roi', 'dmpfc_roi', 'left_vs_roi'};
else
    ROIs = {'vmPFC'};
end

%% Load data
% Load behavioural data
subjList = str2double(participants); 

% Load neural data
timefolder = fullfile(neural_dat);

% Load timings and regressor table
load('mri_timings.mat');
load('regressor_table.mat')

%% 
for iregion = 1:length(ROIs)
    for isub=1:numel(subjList)

        subjInd = find(mri_timings.subjList==subjList(isub));
        onset.subj{isub}.ID = subjList(isub);
        nRuns = numel(mri_timings.onsets.choice_onsets{subjInd,1}); 
        run_data = cell(nRuns,1); % empty cell for saving the epoched file from each run so they can be combined in a single subject 

        out.choiceOnset.Region(iregion).trial_indices{isub}.run{1} = mri_timings.trial_indeces{subjInd}{1}; 
        out.choiceOnset.Region(iregion).trial_indices{isub}.run{2} = mri_timings.trial_indeces{subjInd}{2};
        out.choiceOnset.Region(iregion).subjList = mri_timings.subjList;

        if subjInd == 23
            out.choiceOnset.Region(iregion).trial_indices{isub}.run{3} = mri_timings.trial_indeces{subjInd}{3};
        end

        for irun = 1:nRuns

            if subjInd ~= 23 && irun == 3
                continue
            end
            
            onset.subj{isub}.run{irun} = mri_timings.onsets.choice_onsets{subjInd,1}{irun,1};

            timefile = [timefolder, 'cluster_masks/timeCourse_files/',ROIs{iregion},'/sub_',num2str(subjList(isub), '%02d'),'_run',num2str(irun),'.txt']; 
            timecourseT=load(timefile); % load the timecourse (BOLD responses at the measured times, measurements are one TR apart)
            timecourse = timecourseT - mean(timecourseT); % de-mean 

            
        end
    end
end