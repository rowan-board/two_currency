%% read text files and combine into a regressor table 
% struct containing 
% participants
% runs 
% each regressor

%%
% get relevant paths
behave='C:\Users\Rowan\Desktop\Oxford\behavioural_analysis\';

% set the subjects and the runs
subjects = {'02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23'};
runs = 3; % 3 runs to accomodate patient, just skip unless patient subject number
regressors = 6;

% create empty cell arrays over the number of subjects, runs, and regressors which we will use to make structure 
subs = cell(numel(subjects),1);
r = cell(runs,1);
reg = cell(regressors,1);
x = cell(numel(subjects),1);

% create data_table struct with a cell for each participant
data_table.participants = subs;

% split motivation/devaluation signal into parts?
split = false;

%%
% loop over subjects
for subj = 1:numel(subjects)

    % within each participant cell, place a cell for each run
    data_table.participants{subj}.runs=r;

    % make subject ID string for paths
    subj_id = ['sub_',subjects{subj}];
    
    % loop over runs
    for run = 1:runs
        
        % check if one of the files exists in run 3, if no then skip as no
        % run 3
        if run == 3 && isfile([behave, 'text_files\', subj_id, '\run_3\reward_3.txt'])
            % run 3 exists
        elseif run == 3 && isfile([behave, 'text_files\', subj_id, '\run_3\reward_3.txt']) == false
            % skip run as doesn't exist
            continue
        end
        
        % within each run cell place a cell for each regressor
        data_table.participants{subj}.runs{run}.regressor = reg;

        run_id = ['run_', num2str(run)];
            
        trials = importdata([behave, 'text_files\', subj_id, '\', run_id, '\trial_no_', num2str(run), '.txt']);
        reward = importdata([behave, 'text_files\', subj_id, '\', run_id, '\reward_', num2str(run), '.txt']);
        chosen_ev = importdata([behave, 'text_files\', subj_id, '\', run_id, '\chosen_ev_', num2str(run), '.txt']);
        reward_diff = importdata([behave, 'text_files\', subj_id, '\', run_id, '\reward_diff_', num2str(run), '.txt']);
        chosen_diff = importdata([behave, 'text_files\', subj_id, '\', run_id, '\chosen_diff_', num2str(run), '.txt']);
        motivation = importdata([behave, 'text_files\', subj_id, '\', run_id, '\motivation_', num2str(run), '.txt']);

        data_table.participants{subj}.runs{run}.regressor{1}=trials(:,3);
        data_table.participants{subj}.runs{run}.regressor{2}=reward(:,3);
        data_table.participants{subj}.runs{run}.regressor{3}=chosen_ev(:,3);
        data_table.participants{subj}.runs{run}.regressor{4}=reward_diff(:,3);
        data_table.participants{subj}.runs{run}.regressor{5}=chosen_diff(:,3);
        data_table.participants{subj}.runs{run}.regressor{6}=motivation(:,3);

        if split == true
            m = motivation(:,1);
    
            for i = 1:length(motivation)
                if m(i)==1 && m(i-8)==1
                    m(i)=-1;
                elseif m(i)==1 && m(i-8)==-1
                    m(i)=-1;
                end
            end
    
            extinction=m==1;
            reinstatement=m==-1;
            reinstatement = double(reinstatement);
    
            for k = 1:length(reinstatement)
                if reinstatement(k)==1 && reinstatement(k-7)==1
                    reinstatement(k)=-1;
                elseif reinstatement(k)==1 && reinstatement(k-7)==-1
                    reinstatement(k)=-1;
                end
            end
    
            r1=reinstatement==1;
            r2=reinstatement==-1;
    
            data_table.participants{subj}.runs{run}.regressor{1}=trials(:,3);
            data_table.participants{subj}.runs{run}.regressor{2}=reward(:,3);
            data_table.participants{subj}.runs{run}.regressor{3}=chosen_ev(:,3);
            data_table.participants{subj}.runs{run}.regressor{4}=reward_diff(:,3);
            data_table.participants{subj}.runs{run}.regressor{5}=chosen_diff(:,3);
            data_table.participants{subj}.runs{run}.regressor{6}=extinction;
            data_table.participants{subj}.runs{run}.regressor{7}=r1;
            data_table.participants{subj}.runs{run}.regressor{8}=r2;
    
            x{subj}=numel(trials(:,3));
        end

    end

end

save('regressor_table.mat', 'data_table')