function [out] = upsample_tc_all(timefolder,mri_timings,ROIs,subjList)
% Function to cut the time courses to choice onset 

% Input
    % timefolder: (string) path to folder containing extracted timecourses for ROIs
    % mri_timings: (structure) data structure of onset times for main events 
    % ROIs: (structure) list of strings for ROIs to upsample
    % subjList: array of subject IDs
    
% Output
    % structure with time courses upsampled to choice onset 

% Presettings 

TR=1.23;   % from scan settings
upsample=10; % upsampling factor (interpolate activity in units of 1/upsampling factor) 
norm_timec=1; % normalise time course
    
% Load, cut (according to timings as for FSL) and upsample the time course file

for iregion = 1:length(ROIs)
    
    for isub= 1:numel(subjList) 
        
         
         subjInd = find(mri_timings.subjList==subjList(isub));
         onset.subj{isub}.ID = subjList(isub); 
         
         nRuns = numel(mri_timings.onsets.choice_onsets{subjInd,1}); 
         run_data = cell(nRuns,1); % empty cell for saving the epoched file from each run so they can be combined in a single subject 
         out.choiceOnset.Region(iregion).trial_indices{isub}.run{1} = mri_timings.trial_indeces{subjInd}{1}; 
         out.choiceOnset.Region(iregion).trial_indices{isub}.run{2} = mri_timings.trial_indeces{subjInd}{2};
         out.choiceOnset.Region(iregion).subjList = mri_timings.subjList;
         
        for irun = 1:nRuns
        
            
            onset.subj{isub}.run{irun} = mri_timings.onsets.choice_onsets{subjInd,1}{irun,1};
        
            
            % Load and upsample timecourse
            timefile = [timefolder, '/',ROIs{iregion}, '/timeCourse_files/sub_',num2str(subjList(isub), '%02d'),'_run',num2str(irun),'_all.txt']; 
            timecourseT=load(timefile); % load the timecourse (BOLD responses at the measured times, measurements are one TR apart)
            % crop out first 3 rows as these contain voxel coords
            timecourseT(1:3,:) = [];
            
            % get the number of voxels and the number of measurements taken
            nVoxels=size(timecourseT);
            nTR=nVoxels(1);
            nVoxels=nVoxels(2);

            % pre assign array
            timecourse=zeros(nTR,nVoxels);

            % loop over voxels 
            for i = 1:nVoxels %no. of voxels
                timecourse(:,i) = timecourseT(:,i) - mean(timecourseT(:,i)); % de-mean 
            end

            if norm_timec %normalise the timecourse
                %timecourse=zscore(timecourse);
                for i = 1:nVoxels %no. of voxels
                    timecourse(:,i) = normalize(timecourse(:,i)); 
                end
            end
            
            % pre assign arrays
            ts_ups=cell(1,nVoxels);
            
            % Upsample the time course
            for i = 1:nVoxels %no. of voxels
                t=0:length(timecourse(:,1))-1; % number of time points across run 
                %t_ups=0:(1/upsample):length(timecourse)-1; % number of upsampled time points (interpolating by 1/upsample) 
                t_ups=0:(1/upsample):length(timecourse(:,1)); % upsampled time points (interpolating by 1/upsample) 
                ts_ups{1,i}=spline(t,timecourse(:,i),t_ups); % upsample the timecourse (i.e. the bold values) – interpolate values 
            end
            %ts_upsDerivative=diff([0 ts_ups]);
            
            % Epoch the time course by choice onset 

            window = 12; % epoch window duration in seconds
            TCShift=0; % option for time shift before or after event 
            pseudosamples=round((window./TR)*upsample); % number of pseudosamples in epoch
            
            onsets = round(((onset.subj{isub}.run{irun}+TCShift)./TR)*upsample); % find index of onset within timecourse course data by dividing by TR and multiplying by upsample
            
            tmpmat=repmat(onsets,1,pseudosamples)+repmat(0:pseudosamples-1,size(onsets,1),1); % extend to index timings beyond onsets 
            
            
            % If MRI volumes are missing for the lag period on the final trial, delete the final trial
            while sum(tmpmat>numel(ts_ups{1,1}),'all')>0 % check for any time points which extend beyond MRI time course data 
                tmpmat(end,:) = []; 
                out.choiceOnset.Region(iregion).trial_indices{isub}.run{irun}(end) = []; 
            end
            
            
            %psDerivative=ts_upsDerivative(tmpmat);
            
            % save time course for this run
            run_data{irun}=cell(1,nVoxels);
            for i=1:nVoxels
                run_data{irun}{1,i} = ts_ups{1,i}(tmpmat); % index the upsampled time course data using the onset timings matrix 
            end
            
        end
            out.choiceOnset.Region(iregion).epoched_timecourses{isub}.run = run_data;  
            
            % can't concatenate across runs as different number of voxels
            % across runs

            %out.choiceOnset.Region(iregion).epoched_timecourses{isub}.all_runs = [out.choiceOnset.Region(iregion).epoched_timecourses{isub}.run{1}; out.choiceOnset.Region(iregion).epoched_timecourses{isub}.run{2}];
            %out.choiceOnset.Region(iregion).trial_indices{isub}.all_runs = [out.choiceOnset.Region(iregion).trial_indices{isub}.run{1}, out.choiceOnset.Region(iregion).trial_indices{isub}.run{2}];
            
            out.choiceOnset.Region(iregion).name    = ROIs{iregion};
            out.choiceOnset.info                    = 'Time course time locked to the choice onset phase';
            out.choiceOnset.pseudosamples           = pseudosamples;
            out.choiceOnset.window                  = window;
            out.choiceOnset.TCShift                 = TCShift;
            out.choiceOnset.TR                      = TR;
            out.choiceOnset.upsample                = upsample;



       
    end

end