function [out] = upsample_tc(timefolder,mri_timings,ROIs,subjList)
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

         if subjInd == 23
             out.choiceOnset.Region(iregion).trial_indices{isub}.run{3} = mri_timings.trial_indeces{subjInd}{3};
         end

         
         for irun = 1:nRuns

            if subjInd ~= 23 && irun == 3
                continue
            end
            
            onset.subj{isub}.run{irun} = mri_timings.onsets.choice_onsets{subjInd,1}{irun,1};
        
            
            % Load and upsample timecourse
            timefile = [timefolder, 'cluster_masks/timeCourse_files/', ROIs{iregion}, '/sub_', num2str(subjList(isub), '%02d') '_run', num2str(irun), '.txt'];
            %timefile = [timefolder, ROIs{iregion},'/timeCourse_files/sub_',num2str(subjList(isub), '%02d'),'_run',num2str(irun),'.txt']; 
            timecourseT=load(timefile); % load the timecourse (BOLD responses at the measured times, measurements are one TR apart)
            timecourse = timecourseT - mean(timecourseT); % de-mean 

            if norm_timec %normalise the timecourse
                %timecourse=zscore(timecourse);
                timecourse = normalize(timecourse); 
            end

            % Upsample the time course
            
            t=0:length(timecourse)-1; % number of time points across run 
            %t_ups=0:(1/upsample):length(timecourse)-1; % number of upsampled time points (interpolating by 1/upsample) 
            t_ups=0:(1/upsample):length(timecourse); % upsampled time points (interpolating by 1/upsample) 
            ts_ups=spline(t,timecourse,t_ups); % upsample the timecourse (i.e. the bold values) – interpolate values 
            
            %ts_upsDerivative=diff([0 ts_ups]);
            
            % Epoch the time course by choice onset 

            window = 12; % epoch window duration in seconds
            TCShift=0; % option for time shift before or after event 
            pseudosamples=round((window./TR)*upsample); % number of pseudosamples in epoch
            
            onsets = round(((onset.subj{isub}.run{irun}+TCShift)./TR)*upsample); % find index of onset within timecourse course data by dividing by TR and multiplying by upsample
            
            tmpmat=repmat(onsets,1,pseudosamples)+repmat(0:pseudosamples-1,size(onsets,1),1); % extend to index timings beyond onsets 
            
            
            % If MRI volumes are missing for the lag period on the final trial, delete the final trial
            while sum(tmpmat>numel(ts_ups),'all')>0 % check for any time points which extend beyond MRI time course data 
                tmpmat(end,:) = []; 
                out.choiceOnset.Region(iregion).trial_indices{isub}.run{irun}(end) = []; 
            end
            
            
            %psDerivative=ts_upsDerivative(tmpmat);
            
            % save time course for this run
            run_data{irun} = ts_ups(tmpmat); % index the upsampled time course data using the onset timings matrix 
            
            
         end

        out.choiceOnset.Region(iregion).epoched_timecourses{isub}.run = run_data;  
        
        % concatenate across runs 
        if isub ~= 22
            out.choiceOnset.Region(iregion).epoched_timecourses{isub}.all_runs = [out.choiceOnset.Region(iregion).epoched_timecourses{isub}.run{1}; out.choiceOnset.Region(iregion).epoched_timecourses{isub}.run{2}];
            out.choiceOnset.Region(iregion).trial_indices{isub}.all_runs = [out.choiceOnset.Region(iregion).trial_indices{isub}.run{1}, out.choiceOnset.Region(iregion).trial_indices{isub}.run{2}];
        end
        
        if isub == 22
            out.choiceOnset.Region(iregion).epoched_timecourses{isub}.all_runs = [out.choiceOnset.Region(iregion).epoched_timecourses{isub}.run{1}; out.choiceOnset.Region(iregion).epoched_timecourses{isub}.run{2};  out.choiceOnset.Region(iregion).epoched_timecourses{isub}.run{3}];
            out.choiceOnset.Region(iregion).trial_indices{isub}.all_runs = [out.choiceOnset.Region(iregion).trial_indices{isub}.run{1}, out.choiceOnset.Region(iregion).trial_indices{isub}.run{2}, out.choiceOnset.Region(iregion).trial_indices{isub}.run{3}];
        end

        out.choiceOnset.Region(iregion).name    = ROIs{iregion};
        out.choiceOnset.info                    = 'Time course time locked to the choice onset phase';
        out.choiceOnset.pseudosamples           = pseudosamples;
        out.choiceOnset.window                  = window;
        out.choiceOnset.TCShift                 = TCShift;
        out.choiceOnset.TR                      = TR;
        out.choiceOnset.upsample                = upsample;

       
    end

end