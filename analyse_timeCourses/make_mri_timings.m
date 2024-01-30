%% make mri_timings structure 
% take all of the timings from the mri experiment and turn them into a
% structure that can be used for the time course analysis scripts 

% get relevant paths
behave='C:\Users\Rowan\Desktop\Oxford\behavioural_analysis\';

% structure name is mri_timings
mri_timings.subjList = [01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23]; % list of subject Id's
n = numel(mri_timings.subjList);
runs = 3; 

% lets make an empty structure of subjects and runs then fill with the
% choice data
subs = cell(n,1);
r = cell(runs,1);

subT = cell(n,1);
for i=1:n
    subs{i}=r;
    
    subject = mri_timings.subjList(i);
    subj=num2str(subject, '%02d');

    x = importdata([behave 'text_files/sub_' subj '/run_1/trial_no_1.txt']);
    x = x(:,1);

    y = importdata([behave 'text_files/sub_' subj '/run_2/trial_no_2.txt']);
    y = y(:,1);

    subs{i}{1} = x;
    subs{i}{2} = y;

    subT{i}{1} = 1:numel(x);
    subT{i}{2} = 1:numel(y);

    if isfile([behave 'text_files/sub_' subj '/run_3/trial_no_3.txt'])
        z = importdata([behave 'text_files/sub_' subj '/run_3/trial_no_3.txt']);
        z = z(:,1);
        subs{i}{3}=z;
        subT{i}{3}=1:numel(z);
    end

end
mri_timings.onsets.choice_onsets = subs;
mri_timings.trial_indeces = subT;

save('mri_timings.mat', 'mri_timings')
