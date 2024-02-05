%% Load SNIRF pipeline processed data

load('../Data/mat_files/SubjStats_setup_1_cbsi_vardpf.mat');
SubjStats = SubjStats_cbsi;
demographics = nirs.createDemographicsTable(SubjStats);
disp(demographics);

%% Set groups in demographics

demographics = nirs.createDemographicsTable(SubjStats);
disp(demographics);

ya_ids = readtable("../Data/identifiers_YA.csv");
oa_ids = readtable("../Data/identifiers_OA.csv");
pd_ids = readtable("../Data/identifiers_PD.csv");
subj_ids = cell2table(demographics.SubjectID, "VariableNames", ["id_nummer"]);

ya_idx = ismember(subj_ids, ya_ids);
oa_idx = ismember(subj_ids, oa_ids);
pd_idx = ismember(subj_ids, pd_ids);

demographics.group = repmat("NA", height(demographics), 1);
demographics.group(ya_idx) = "YA";
demographics.group(oa_idx) = "OA";
demographics.group(pd_idx) = "PD";

job=nirs.modules.AddDemographics;
job.varToMatch = 'UUID';
job.demoTable=demographics;
SubjStats=job.run(SubjStats);

demographics = nirs.createDemographicsTable(SubjStats);
disp(demographics);

fprintf('Sum YA: %d\n', sum(ya_idx));
fprintf('Sum OA: %d\n', sum(oa_idx));
fprintf('Sum PD: %d\n', sum(pd_idx));

%% Add covariates

% Get the demographics
balance_data = readtable('../Data/REDcap_data/MiniBEST_data.csv'); 
w12_data = readtable('../Data/REDcap_data/Walk12_data.csv'); 
hads_data = readtable('../Data/REDcap_data/HADS_data.csv'); 
gait_data = readtable('../Data/temp_data/validation_data_imu_gait_parameters.csv'); 
updrs_data = readtable('../Data/REDcap_data/UPDRS_data.csv'); 

% Prepare the table
demographics.hy = repmat("NA", height(demographics), 1);
demographics.updrs_3_motor = NaN(height(demographics),1);
demographics.balance = NaN(height(demographics),1);
demographics.w12 = NaN(height(demographics),1);
demographics.hads_anxiety = NaN(height(demographics),1);
demographics.st_walk_speed = NaN(height(demographics),1);
demographics.dt_cost_walk_speed = NaN(height(demographics),1);
demographics.st_step_time_var = NaN(height(demographics),1);

% Add balance data
for idx=1:height(balance_data)
    subj_id_seek = string(balance_data.id_nummer(idx)); 
    match_idx = strcmp(demographics.SubjectID, subj_id_seek);
    if sum(match_idx) < 1
        continue
    end
    demographics(match_idx,:).balance = balance_data(idx,:).mb_total;
end

% W12
for idx=1:height(w12_data)
    subj_id_seek = string(w12_data.id_nummer(idx)); 
    match_idx = strcmp(demographics.SubjectID, subj_id_seek);
    if sum(match_idx) < 1
        continue
    end
    demographics(match_idx,:).w12 = w12_data(idx,:).g12_sum;
end

% HADS
for idx=1:height(hads_data)
    subj_id_seek = string(hads_data.id_nummer(idx)); 
    match_idx = strcmp(demographics.SubjectID, subj_id_seek);
    if sum(match_idx) < 1
        continue
    end
    demographics(match_idx,:).hads_anxiety = hads_data(idx,:).sum_anxiety;
end

% Gait data
for idx=1:height(gait_data)
    subj_id_seek = string(gait_data.subject(idx)); 
    match_idx = strcmp(demographics.SubjectID, subj_id_seek);
    if sum(match_idx) < 1
        continue
    end
    demographics(match_idx,:).st_walk_speed = gait_data(idx,:).st_walk_speed;
    demographics(match_idx,:).dt_cost_walk_speed = gait_data(idx,:).dt_cost_walk_speed;
    demographics(match_idx,:).st_step_time_var = gait_data(idx,:).st_step_time_var;
end

% UPDRS3 motor score and HY
for idx=1:height(updrs_data)
    subj_id_seek = string(updrs_data.id_nummer(idx)); 
    match_idx = strcmp(demographics.SubjectID, subj_id_seek);
    if sum(match_idx) < 1
        continue
    end
    updrs_3_scores = table2array(updrs_data(idx,33:65));
    if sum(isnan(updrs_3_scores)) == length(updrs_3_scores)
        continue
    end
    updrs_3 = sum(updrs_3_scores, "omitnan");
    demographics(match_idx,:).updrs_3_motor = updrs_3;
    hy_value = updrs_data(idx,:).mdsupdrs3_hy;
    demographics(match_idx,:).hy(hy_value == 1 | hy_value == 2) = "HY_1_2";
    demographics(match_idx,:).hy(hy_value == 3 | hy_value == 4) = "HY_3_4";
end

job = nirs.modules.AddDemographics;
job.demoTable = demographics;
job.varToMatch = 'UUID';
SubjStats = job.run(SubjStats);

%% ROI 

% Set up ROI
source = [NaN NaN NaN NaN NaN NaN NaN]';    
detector = [1 2 3 4 5 6 7]';
ROI_PFC = table(source,detector);

% ba9 and ba46
source = [4 7 4 1 4 6 8 3 2]';
detector = [5 5 2 2 4 7 7 1 1]';
ROI_ba_9_46 = table(source,detector);
weights_ba9 = [0.689, 0.684, 0.689, 0.684, 0.618]';
weights_ba9 = normalize(weights_ba9, 'norm', 1);
weights_ba46 = [0.47, 0.43, 0.47, 0.43]';
weights_ba46 = normalize(weights_ba46, 'norm', 1);
weights = [weights_ba9; weights_ba46];
weights = normalize(weights, 'norm', 1);
ROI_ba_9_46.weight = weights;

%% Write betas

ROIs = {ROI_ba_9_46};
ROI_names = ["BA9+46"];

for i=1:length(ROIs)
    
    % Get ROI results for each condition
    roi_result_all = nirs.util.roiAverage(SubjStats, ROIs{i}, ROI_names(i));
    roi_result_all.SubjectID = repmat('FNPXXXX',height(roi_result_all),1);
    for idx=1:height(roi_result_all)
        subj_idx = str2num(roi_result_all.FileIdx{idx});
        roi_result_all.SubjectID(idx,:) = SubjStats(subj_idx).demographics.SubjectID;
        roi_result_all.group(idx,:) = SubjStats(subj_idx).demographics.group;
        roi_result_all.st_walk_speed(idx,:) = SubjStats(subj_idx).demographics.st_walk_speed;
        roi_result_all.st_step_time_var(idx,:) = SubjStats(subj_idx).demographics.st_step_time_var;
    end
    fname = 'subject_level_ROI_beta_cbsi_' + ROI_names(i) + '_ALL.csv';
    writetable(roi_result_all, fname);
    
end
