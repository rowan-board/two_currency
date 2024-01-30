function [betas,SE, tStat] = fit_tc_GLM(reg_mat, timedat, trials_idx, intercept, nRegressors_int) 

ntimep = size(timedat,2); 

betas = nan(ntimep,nRegressors_int); %subs, region, timepoints, regressors
tStat = nan(ntimep,nRegressors_int); %subs, region, timepoints, regressors
SE = nan(ntimep,nRegressors_int); %subs, region, timepoints, regressors

 for tp = 1:ntimep
    
    % Index by trials and time point 
	ydat = timedat(trials_idx,tp);
	pred = reg_mat(trials_idx,:);
    
    % Fit to each time point for that participant 

            if intercept
                modelFit = fitglm(pred,ydat);
                betas(tp,:) = modelFit.Coefficients.Estimate(1:end); 
                SE(tp,:) = modelFit.Coefficients.SE(1:end); 
                tStat(tp,:) = modelFit.Coefficients.tStat(1:end); 
            else
                modelFit = fitglm(pred,ydat,'Intercept',false); 
                betas(tp,:) = modelFit.Coefficients.Estimate(1:end); 
                SE(tp,:) = modelFit.Coefficients.SE(1:end); 
                tStat(tp,:) = modelFit.Coefficients.tStat(1:end); 
            end
            
 end    