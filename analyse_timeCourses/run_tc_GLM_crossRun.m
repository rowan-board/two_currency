function [betas,SE, tStat] = time_course_regression(timec, subjList, aggregateTable,  regressors,intercept, treat_regressors, analysis_regions, trial_condition_agg)
% runs for each subject, for each time point a separate regression with all
% regressors 
% treat_regressors – 0=unnormalised values, 1=normalise, 2 = demean (but
% don't divide by variance) 

% If no index supplied, use full table
if ~exist('trial_condition_agg','var')
    trial_condition_agg = logical(ones(height(aggregateTable),1)); 
end
if ~exist('intercept','var')
    intercept = 0;
end

% Define number of regressors and time points 
nRegressors = numel(regressors);
if intercept 
    nRegressors_int = nRegressors + 1;
else
    nRegressors_int = nRegressors;
end
nregions = numel(analysis_regions); 
ntimep = size(timec.Region(1).epoched_timecourses{1}.all_runs,2);

 
% Preassign arrays 
betas = nan(length(subjList),nregions, ntimep,nRegressors_int); %subs, region, timepoints, regressors
tStat = nan(length(subjList),nregions, ntimep,nRegressors_int); %subs, region, timepoints, regressors
SE = nan(length(subjList),nregions, ntimep,nRegressors_int); %subs, region, timepoints, regressors

for reg = 1:nregions
    

    for p = 1:length(subjList)
        
        %subjInd = find(GLM.subjList==subjList(isub)); 
        
        timedat = timec.Region(reg).epoched_timecourses{p}.all_runs; 
        trialidx = timec.Region(reg).trial_indices{p}.all_runs; 

        % load subject specific data table 
        data_table = aggregateTable(aggregateTable.participant==subjList(p),:); % subject specific data table 
        idx_trials = trial_condition_agg(aggregateTable.participant==subjList(p)); % subject specific trial condition indices 
        idx_trials = idx_trials(trialidx); % use mri trial indices to crop the logical condition indices to fit quantity of mri data 
        nTrials = numel(trialidx); 

        % make regressor table
        regMat = zeros(nTrials,nRegressors); 
        for r = 1:nRegressors
            regs = eval(strcat('data_table.',regressors{r}));
            regMat(:,r) = regs(trialidx); 
            
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
        end


       [betas(p,reg,:,:),SE(p,reg,:,:), tStat(p,reg,:,:)] = fit_tc_GLM(regMat, timedat, idx_trials,intercept,nRegressors_int); 
       
    end

end
