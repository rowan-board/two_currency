function [betas,SE, tStat] = time_course_regression(choiceOnset, subjList, data_table,  regressors, intercept, treat_regressors, analysis_regions)
% runs for each subject, for each time point a separate regression with all
% regressors 
% treat_regressors – 0=unnormalised values, 1=normalise, 2 = demean (but
% don't divide by variance) 

% If no index supplied, use full table
%if ~exist('trial_condition_agg','var')
%    trial_condition_agg = logical(ones(height(aggregateTable),1)); 
%end
%if ~exist('intercept','var')
%    intercept = 0;
%end

% Define number of regressors and time points 
nRegressors = numel(regressors);
if intercept 
    nRegressors_int = nRegressors + 1;
else
    nRegressors_int = nRegressors;
end
nregions = numel(analysis_regions); 
ntimep = size(choiceOnset.Region(1).epoched_timecourses{1}.all_runs,2);

 
% Preassign arrays 
betas = nan(length(subjList),nregions, ntimep,nRegressors_int); %subs, region, timepoints, regressors
tStat = nan(length(subjList),nregions, ntimep,nRegressors_int); %subs, region, timepoints, regressors
SE = nan(length(subjList),nregions, ntimep,nRegressors_int); %subs, region, timepoints, regressors

for reg = 1:nregions
    
    for p = 1:length(subjList)
    
        
        % check if third run exists, if no just run 2 runs (mainly a check
        % for patient as they have 3 runs)
        if choiceOnset.Region(reg).epoched_timecourses{p}.run{3}
            nRuns = numel(choiceOnset.Region(reg).epoched_timecourses{p}.run);
        else
            nRuns = 2;
        end

        
        beta_run = NaN(nRuns, ntimep,nRegressors_int); 
        tStat_run = NaN(nRuns, ntimep,nRegressors_int); 
        SE_run = NaN(nRuns, ntimep,nRegressors_int); 
        
        for irun = 1:nRuns        
        
            timedat = choiceOnset.Region(reg).epoched_timecourses{p}.run{irun}; 
            trialidx = choiceOnset.Region(reg).trial_indices{p}.run{irun};

            % load subject specific data table
            if irun==1
              x = (p*2)-1;
            else
              x = (p*2);
            end

            regressor_table = data_table.participants{p}.runs{irun}.regressor(:,:); % subject and run specific data table 
            idx_trials = 1:numel(regressor_table{1}); % subject and run specific trial condition indices 
            idx_trials = idx_trials(trialidx); % use mri trial indices to crop the logical condition indices to fit quantity of mri data 
            nTrials = numel(trialidx); 

            % edit regressor table so that each trial has it's own cell,
            % instead of each regressor has 1 cell
            reg_mat = zeros(nTrials,nRegressors);

            % If MRI volumes are missing for the lag period on the final
            % trial, delete final trials
            % check for any time points which extend beyond MRI time course data 
            x=numel(reg_mat(:,1));

            for z = 1:nRegressors
                reg_mat(:,z)=regressor_table{z}(1:x);
            end


            % regressor treatment, will need editing if you want to use
            if treat_regressors ==1
                regMat(:,r) = normalize(regMat(:,r)); 
            elseif treat_regressors==2
                regMat(:,r) = regMat(:,r) - mean(regMat(:,r)); 
            elseif treat_regressors==3
                nonzeros = regMat(:,r)~=0;
                regMat(nonzeros,r) = normalize(regMat(nonzeros,r));
            elseif treat_regressors==4
                regMat(:,r) = normalize(regMat(:,r));
                regMat(isnan(regMat(:,r)),r) = 0;
            end


           [beta_run(irun, :,:), SE_run(irun, :,:), tStat_run(irun, :,:)] = fit_tc_GLM(reg_mat, timedat, idx_trials,intercept,nRegressors_int); 
           
        end
        
        betas(p,reg,:,:) = mean(beta_run,1); 
        SE(p,reg,:,:)= mean(SE_run,1); 
        tStat(p,reg,:,:)= mean(tStat_run,1); 
    end

end
