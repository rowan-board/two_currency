%% Define paths 

homepath = '/Users/Rowan/Desktop/oxford';
addpath(genpath(homepath)); 

neural_dat = [homepath '/analyse_timeCourses/'];

% voxel wise (1) or mean ROI (0)?
voxel = 0;

participants = {'02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23'}; 

if voxel==0
    ROIs = {'vmPFC','dmPFC','VS','mFP','lFP'};
    ROIs = {'right_vs_roi', 'left_vs_roi', 'ventral_striatum_roi', 'left_amygdala_roi'};
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
load('regressor_table.mat');


%% Upsample time course by event

if voxel == 0
    upsampledData = upsample_tc(timefolder,mri_timings,ROIs,subjList); 
else
    upsampledData = upsample_tc_all(timefolder,mri_timings,ROIs,subjList);
end
choiceOnset = upsampledData.choiceOnset;
%% Run the regression
regressors= {'trials','reward','chosen_ev','reward_diff','chosen_diff','motivation'};
%regressors= {'trials','reward','chosen_ev','reward_diff','chosen_diff','extinction','r1','r2'};
% Run model

nRegressors = numel(regressors); 
treat_regressors=0; %  0=unnormalised values, 1=normalise, 2 = demean (but
% don't divide by variance), 3=normalise but ignoring all zeros 
analysis_regions = 1:numel(ROIs); 


nModels = 1;
models = cell(nModels,1);
intercept = 0;

if voxel == 0
    [models{1}.betas,models{1}.SE, models{1}.tStat] = run_tc_GLM_byRun(choiceOnset, subjList, data_table,  regressors,intercept, treat_regressors, analysis_regions);
    models{1}.name = 'all trials';
else
    [models{1}.betas,models{1}.SE, models{1}.tStat] = run_tc_GLM_byVoxel(choiceOnset, subjList, data_table,  regressors,intercept, treat_regressors, analysis_regions);
    models{1}.name = 'all trials';
end

%% plot regression
nModels = numel(models); 
ntimep = size(choiceOnset.Region(1).epoched_timecourses{1}.all_runs,2);   % number of timepoints in each epoch
modelColours = {'r','b','g','k'};



% regressorNames = {'Trials','Objective Reward','Objective EV','Devaluation of Reward','Devaluation of EV','Extinction','Reinstatement (1st Half)','Reinstatement (2nd Half)'};
regressorNames = {'Trials','Objective Reward','Objective EV','Devaluation of Reward','Devaluation of EV','Motivation'};

for r = 1:numel(regressorNames)
    regressorNames{r} = strrep(regressorNames{r},'_',' ');
end

regionNames = [];
for reg = 1:numel(analysis_regions)
    regionNames = [regionNames, {strrep(choiceOnset.Region(reg).name,'_',' ')}];
end

changeColour = 1;
for reg = 1:numel(analysis_regions)
    
    figure; 
    sgtitle(regionNames(reg)); 
    
    count = 0;
    for b = 1:nRegressors
        
        count = count+1;
        subplot(4,4,count);
        title(regressorNames{b});  

        for m = 1:nModels

           beta_weights = models{1}.betas; 
           standarderror = std(beta_weights(:,reg,:,b))./sqrt(length(subjList)); 
           y = mean(beta_weights(:,reg,:,b)); 
           x = [choiceOnset.window/ntimep:choiceOnset.window/ntimep:choiceOnset.window]; 
           xlabel('time (s)');
           ylabel('beta coefficient');
           
           shadedErrorBar(x,y,standarderror,'lineProps',{modelColours{1}});

           % if devaluation of reward in vmPFC then overlay devaluation of
           % chosen
           if reg == 1 & b == 5
               hold on 
               standarderror = std(beta_weights(:,reg,:,b))./sqrt(length(subjList)); 
               y = mean(beta_weights(:,reg,:,b-2));
               shadedErrorBar(x,y,standarderror,'lineProps',{modelColours{2}});
               legend({'DevEV','EV'},'location','northeastoutside','Orientation','vertical','AutoUpdate','off')
               hold off
           end

           % if devaluation of reward in vS then overlay on reward
           if reg == 3 & b == 4
               hold on 
               standarderror = std(beta_weights(:,reg,:,b))./sqrt(length(subjList)); 
               y = mean(beta_weights(:,reg,:,b-2));
               shadedErrorBar(x,y,standarderror,'lineProps',{modelColours{2}});
               legend({'DevRwd','Rwd'},'location','northeastoutside','Orientation','vertical','AutoUpdate','off')
               hold off
           end

           xticks([0:2:choiceOnset.window])
           xlim([0, 12]); 
           hold on; 
        end

        % set title as axis labels

        if reg == 3 & b == 4
            title('Reward and Devaluation of Reward'); 
        elseif reg == 1 & b == 5
            title('Objective EV and EV Devaluation'); 
        else
            title(regressorNames{b});
        end
        yline(0,'k'); 
        xlabel('time (s)');
        ylabel('mean signal change (a.u)');

    end
end

%% plot regression of patient vs controls
nModels = numel(models); 
ntimep = size(choiceOnset.Region(1).epoched_timecourses{1}.all_runs,2);  % number of timepoints in each epoch
modelColours = {'r','b','g','k'};

regressorNames = {'Trials','Objective Reward','Objective EV','Devaluation of Reward','Devaluation of EV','Motivation'};

for r = 1:numel(regressorNames)
    regressorNames{r} = strrep(regressorNames{r},'_',' ');
end

regionNames = [];
for reg = 1:numel(analysis_regions)
    regionNames = [regionNames, {strrep(choiceOnset.Region(reg).name,'_',' ')}];
end

changeColour = 1;
for reg = 1:numel(analysis_regions)
    
    figure; 
    sgtitle(regionNames(reg)); 
    
    count = 0;
    for b = 1:nRegressors
        
        count = count+1;
        subplot(4,4,count);
        title(regressorNames{b});  

        for m = 1:nModels

           beta_weights = models{1}.betas; 
           standarderror = std(beta_weights(1:21,reg,:,b))./sqrt(length(subjList)); 
           y = mean(beta_weights(1:21,reg,:,b)); 
           x = [choiceOnset.window/ntimep:choiceOnset.window/ntimep:choiceOnset.window]; 
           xlabel('time (s)');
           ylabel('beta coefficient');
           
           shadedErrorBar(x,y,standarderror,'lineProps',{modelColours{1}});

           %vs, devaluation of reward, patient vs control
           if reg == 2 & b == 2
               hold on 
               standarderror = beta_weights(22,reg,:,b)./sqrt(length(subjList)); 
               y = beta_weights(22,reg,:,b);
               shadedErrorBar(x,y,standarderror,'lineProps',{modelColours{2}});
               legend({'ctrl','pat'},'location','northeastoutside','Orientation','vertical','AutoUpdate','off')
               hold off
           end

           %vs, devaluation of reward, patient vs control
           if reg == 2 & b == 4
               hold on 
               standarderror = beta_weights(22,reg,:,b)./sqrt(length(subjList)); 
               y = beta_weights(22,reg,:,b);
               shadedErrorBar(x,y,standarderror,'lineProps',{modelColours{2}});
               legend({'ctrl','pat'},'location','northeastoutside','Orientation','vertical','AutoUpdate','off')
               hold off
           end

           %vs, devaluation of ev, patient vs control
           if reg == 2 & b == 5
               hold on 
               standarderror = beta_weights(22,reg,:,b)./sqrt(length(subjList)); 
               y = beta_weights(22,reg,:,b);
               shadedErrorBar(x,y,standarderror,'lineProps',{modelColours{2}});
               legend({'ctrl','pat'},'location','northeastoutside','Orientation','vertical','AutoUpdate','off')
               hold off
           end

           %vmpfc, devaluation of rwd, patient vs control
           if reg == 1 & b == 4
               hold on 
               standarderror = beta_weights(22,reg,:,b)./sqrt(length(subjList)); 
               y = beta_weights(22,reg,:,b);
               shadedErrorBar(x,y,standarderror,'lineProps',{modelColours{2}});
               legend({'ctrl','pat'},'location','northeastoutside','Orientation','vertical','AutoUpdate','off')
               hold off
           end

           %vmpfc, devaluation of ev, patient vs control
           if reg == 1 & b == 5
               hold on 
               standarderror = beta_weights(22,reg,:,b)./sqrt(length(subjList)); 
               y = beta_weights(22,reg,:,b);
               shadedErrorBar(x,y,standarderror,'lineProps',{modelColours{2}});
               legend({'ctrl','pat'},'location','northeastoutside','Orientation','vertical','AutoUpdate','off')
               hold off
           end

           %vmpfc, m, patient vs control
           if reg == 1 & b == 6
               hold on 
               standarderror = beta_weights(22,reg,:,b)./sqrt(length(subjList)); 
               y = beta_weights(22,reg,:,b);
               shadedErrorBar(x,y,standarderror,'lineProps',{modelColours{2}});
               legend({'ctrl','pat'},'location','northeastoutside','Orientation','vertical','AutoUpdate','off')
               hold off
           end

           %dmpfc, devaluation of ev, patient vs control
           if reg == 3 & b == 5
               hold on 
               standarderror = beta_weights(22,reg,:,b)./sqrt(length(subjList)); 
               y = beta_weights(22,reg,:,b);
               shadedErrorBar(x,y,standarderror,'lineProps',{modelColours{2}});
               legend({'ctrl','pat'},'location','northeastoutside','Orientation','vertical','AutoUpdate','off')
               hold off
           end

           %left vs, devaluation of ev, patient vs control
           if reg == 4 & b == 5
               hold on 
               standarderror = beta_weights(22,reg,:,b)./sqrt(length(subjList)); 
               y = beta_weights(22,reg,:,b);
               shadedErrorBar(x,y,standarderror,'lineProps',{modelColours{2}});
               legend({'ctrl','pat'},'location','northeastoutside','Orientation','vertical','AutoUpdate','off')
               hold off
           end

           %left vs, reward, patient vs control
           if reg == 4 & b == 2
               hold on 
               standarderror = beta_weights(22,reg,:,b)./sqrt(length(subjList)); 
               y = beta_weights(22,reg,:,b);
               shadedErrorBar(x,y,standarderror,'lineProps',{modelColours{2}});
               legend({'ctrl','pat'},'location','northeastoutside','Orientation','vertical','AutoUpdate','off')
               hold off
           end

           %left vs, dev reward, patient vs control
           if reg == 4 & b == 4
               hold on 
               standarderror = beta_weights(22,reg,:,b)./sqrt(length(subjList)); 
               y = beta_weights(22,reg,:,b);
               shadedErrorBar(x,y,standarderror,'lineProps',{modelColours{2}});
               legend({'ctrl','pat'},'location','northeastoutside','Orientation','vertical','AutoUpdate','off')
               hold off
           end


           xticks([0:2:choiceOnset.window])
           xlim([0, 12]); 
           hold on; 
        end

        % set title as axis labels
        title(regressorNames{b});
        yline(0,'k'); 
        xlabel('time (s)');
        ylabel('mean signal change (a.u)');

    end
end



%% extract individual beta weights and plot averages
%regressorNames_nT={'Objective Reward','Objective chosen EV','Devaluation of Reward','Devaluation of EV','Extinction','Reinstatement (1st Half)','Reinstatement (2nd Half)'};
regressorNames_nT = {'Objective Reward','Objective EV','Devaluation of Reward','Devaluation of EV','Motivation'};
ROIs_nT = {'Medial Prefrontal Cortex','Dorsomedial Prefrontal Cortex','Ventral Striatum','Medial Frontal Pole','Lateral Frontal Pole'};


subj_betas = zeros(length(subjList),nRegressors,length(regionNames));



for i=1:length(subjList)
    for j=1:nRegressors
        for k=1:length(regionNames)
            subj_betas(i,j,k)=mean(beta_weights(i,k,:,j));
        end
    end
end

for k=1:length(regionNames)
    
    figure;
    b=bar(mean(subj_betas(:,2:6,k)));
    hold on
    er=errorbar(mean(subj_betas(:,2:6,k)),(std(subj_betas(:,2:6,k))/sqrt(length(subj_betas(:,2:6,k)))));
    er.Color = [0 0 0];                            
    er.LineStyle = 'none'; 
    hold off

    title([ROIs_nT{k}],'FontSize',14); 
    yline(0,'k'); 
    xlabel('regressor','FontSize',20);
    ylabel('average beta coefficient','FontSize',20);
    xticklabels(regressorNames_nT)
    ax = gca;
    ax.XAxis.FontSize = 14; 
    ax.YAxis.FontSize = 14; 

end

% now lets save the averages for each individual in csv files

for i=1:nRegressors
    for j=1:length(regionNames)
        tmp=subj_betas(:,i,j);
        writematrix(tmp, strjoin([neural_dat, string(ROIs(j)), '/', string(regressors(i)), '_', string(ROIs(j)), '.csv'],''));
    end
end

%% extract beta weights and plot averages for patient vs control
regressorNames_nT = {'Objective Reward','Objective EV','Devaluation of Reward','Devaluation of EV','Motivation'};
ROIs_nT = {'right vs roi', 'left vs roi', 'ventral striatum roi', 'left amygdala roi'};

subj_betas = zeros(length(subjList),nRegressors,length(regionNames));

beta_weights = models{1}.betas; 

for i=1:length(subjList)
    for j=1:nRegressors
        for k=1:length(regionNames)
            subj_betas(i,j,k)=mean(beta_weights(i,k,:,j));
        end
    end
end

for k=1:length(regionNames)

    figure;
    b=bar(mean(subj_betas(1:21,2:6,k)));
    hold on
    for i=2:nRegressors
        p=plot(i-1,mean(subj_betas(22,i,k)), '.', Color='red', MarkerSize=10);
        %errorbar(i-1,mean(subj_betas(22,i,k)),(mean(models{1}.betas(22,k,:,i))/mean(models{1}.tStat(22,k,:,i))/2), Color='red');
    end
    er=errorbar(mean(subj_betas(1:21,2:6,k)),(std(subj_betas(1:21,2:6,k))/sqrt(length(subj_betas(1:21,2:6,k)))));
    %for i=2:nRegressors
    %    er_c=errorbar(i,mean(subj_betas(22,i,k)),(mean(models{1}.betas(22,k,:,i))/mean(models{1}.tStat(22,k,:,i))/2), Color='black');
    %end
    er.Color = [0 0 0];                            
    er.LineStyle = 'none'; 
    hold off
        
    

    title([ROIs_nT{k}],'FontSize',14); 
    yline(0,'k'); 
    xlabel('regressor','FontSize',20);
    ylabel('average beta coefficient','FontSize',20);
    xticklabels(regressorNames_nT)
    ax = gca;
    ax.XAxis.FontSize = 14; 
    ax.YAxis.FontSize = 14; 

end

for k=1:length(regionNames)
    ctrl=mean(subj_betas(1:21,2:6,k));
    ctrl_e=std(subj_betas(1:21,2:6,k))/sqrt(length(subj_betas(1:21,2:6,k)));
    pat=[0 0 0 0 0];
    pat_e=[0 0 0 0 0];
    for i=2:nRegressors
        pat(i-1)=mean(subj_betas(22,i,k));
        pat_e(i-1)=mean(models{1}.betas(22,k,:,i))/mean(models{1}.tStat(22,k,:,i))/2;
    end

    figure;
    b=bar([1,2,3,4,5],[ctrl;pat]);

    hold on
    er=errorbar([0.86,1.86,2.86,3.86,4.86],mean(subj_betas(1:21,2:6,k)),(std(subj_betas(1:21,2:6,k))/sqrt(length(subj_betas(1:21,2:6,k)))));
    for i=2:nRegressors
        er_c=errorbar(i-0.86,mean(subj_betas(22,i,k)),(mean(models{1}.betas(22,k,:,i))/mean(models{1}.tStat(22,k,:,i))/2), Color='black');
    end
    er.Color = [0 0 0];                            
    er.LineStyle = 'none'; 
    hold off

    title([ROIs_nT{k}],'FontSize',14); 
    yline(0,'k'); 
    xlabel('regressor','FontSize',20);
    ylabel('average beta coefficient','FontSize',20);
    xticklabels(regressorNames_nT)
    ax = gca;
    ax.XAxis.FontSize = 14; 
    ax.YAxis.FontSize = 14; 

    legend(['Control';'Patient'],Location="northeastoutside")
end

% patient vs control csv roi
% now lets save the averages for each individual in csv files

for i=1:nRegressors
    for j=1:length(regionNames)
        tmp=subj_betas(:,i,j);
        writematrix(tmp, strjoin([neural_dat, string(ROIs(j)), '/', string(regressors(i)), '_', string(ROIs(j)), '.csv'],''));
    end
end

