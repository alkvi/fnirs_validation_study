%% Select dataset

clear;

% One of cbsi, hbo, hbr, hbd, hbt
hemo_measure = "cbsi";

% File locations
subjstats_file = "../Data/mat_files/SubjStats_setup_1_" + hemo_measure + "_vardpf.mat";
results_file = "../Data/temp_data/hypothesis_table_" + hemo_measure + ".csv";

% Load
SubjStats = importdata(subjstats_file);

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

% Normalize covariates for comparable betas
demographics(:,37:43) = normalize(demographics(:,37:43), 'norm');

job = nirs.modules.AddDemographics;
job.demoTable = demographics;
job.varToMatch = 'UUID';
SubjStats = job.run(SubjStats);

%% Select groups

demographics = nirs.createDemographicsTable(SubjStats);

ya_idx = strcmpi(demographics.group,'YA');
oa_idx = strcmpi(demographics.group,'OA');
pd_idx = strcmpi(demographics.group,'PD');

% Keep all subjects
nan_idx = logical(zeros(height(demographics),1));

% Only YA
selected_idx = zeros(size(SubjStats,1),1);
selected_idx(ya_idx,1) = 1;
selected_idx(nan_idx,1) = 0;
selected_idx = logical(selected_idx);
SubjStats_ya = SubjStats(selected_idx);
demographics_ya = nirs.createDemographicsTable(SubjStats_ya);

% Only OA
selected_idx = zeros(size(SubjStats,1),1);
selected_idx(oa_idx,1) = 1;
selected_idx(nan_idx,1) = 0;
selected_idx = logical(selected_idx);
SubjStats_oa = SubjStats(selected_idx);
demographics_oa = nirs.createDemographicsTable(SubjStats_oa);

% Only PD
selected_idx = zeros(size(SubjStats,1),1);
selected_idx(pd_idx,1) = 1;
selected_idx(nan_idx,1) = 0;
selected_idx = logical(selected_idx);
SubjStats_pd = SubjStats(selected_idx);
demographics_pd = nirs.createDemographicsTable(SubjStats_pd);

fprintf('Sum YA: %d\n', length(SubjStats_ya));
fprintf('Sum OA: %d\n', length(SubjStats_oa));
fprintf('Sum PD: %d\n', length(SubjStats_pd));

%% Main condition ROI model

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

% Run group model
job = nirs.modules.MixedEffects();
job.formula = 'beta ~ -1 + cond + (1|SubjectID)'; 
job.dummyCoding = 'full';
job.include_diagnostics = true;
GroupStats_ya = job.run(SubjStats_ya);
GroupStats_oa = job.run(SubjStats_oa);
GroupStats_pd = job.run(SubjStats_pd);

% Collect results
% YA
roi_result_ya = nirs.util.roiAverage(GroupStats_ya, ROI_ba_9_46, 'BA9_46');
AIC = roi_result_ya.model{1}.ModelCriterion.AIC;
BIC = roi_result_ya.model{1}.ModelCriterion.BIC;
roi_result_ya.formula = repmat("beta ~ -1 + cond + (1|SubjectID)", size(roi_result_ya,1),1);
roi_result_ya.AIC = repmat(AIC, size(roi_result_ya,1),1);
roi_result_ya.BIC = repmat(BIC, size(roi_result_ya,1),1);
roi_result_ya.group = repmat("YA", size(roi_result_ya,1),1);

% OA
roi_result_oa = nirs.util.roiAverage(GroupStats_oa, ROI_ba_9_46, 'BA9_46');
AIC = roi_result_oa.model{1}.ModelCriterion.AIC;
BIC = roi_result_oa.model{1}.ModelCriterion.BIC;
roi_result_oa.formula = repmat("beta ~ -1 + cond + (1|SubjectID)", size(roi_result_oa,1),1);
roi_result_oa.AIC = repmat(AIC, size(roi_result_oa,1),1);
roi_result_oa.BIC = repmat(BIC, size(roi_result_oa,1),1);
roi_result_oa.group = repmat("OA", size(roi_result_oa,1),1);

% PD
roi_result_pd = nirs.util.roiAverage(GroupStats_pd, ROI_ba_9_46, 'BA9_46');
AIC = roi_result_pd.model{1}.ModelCriterion.AIC;
BIC = roi_result_pd.model{1}.ModelCriterion.BIC;
roi_result_pd.formula = repmat("beta ~ -1 + cond + (1|SubjectID)", size(roi_result_pd,1),1);
roi_result_pd.AIC = repmat(AIC, size(roi_result_pd,1),1);
roi_result_pd.BIC = repmat(BIC, size(roi_result_pd,1),1);
roi_result_pd.group = repmat("PD", size(roi_result_pd,1),1);

% Put in a single table
main_condition_result = [roi_result_ya; roi_result_oa; roi_result_pd];

% Condition contrasts
condition_contrast_table = [];
c = [0 -1 1];

% YA
ContrastStats = GroupStats_ya.ttest(c);
roi_result = nirs.util.roiAverage(ContrastStats, ROI_ba_9_46, 'BA9_46');
AIC = roi_result_ya.model{1}.ModelCriterion.AIC;
BIC = roi_result_ya.model{1}.ModelCriterion.BIC;
roi_result.formula = repmat("beta ~ -1 + cond + (1|SubjectID)", size(roi_result,1),1);
roi_result.AIC = repmat(AIC, size(roi_result,1),1);
roi_result.BIC = repmat(BIC, size(roi_result,1),1);
roi_result.group = repmat("YA", size(roi_result,1),1);
roi_result.model = repmat({roi_result_ya.model{1}}, size(roi_result,1),1);
condition_contrast_table = [condition_contrast_table; roi_result];

% OA
ContrastStats = GroupStats_oa.ttest(c);
roi_result = nirs.util.roiAverage(ContrastStats, ROI_ba_9_46, 'BA9_46');
AIC = roi_result_oa.model{1}.ModelCriterion.AIC;
BIC = roi_result_oa.model{1}.ModelCriterion.BIC;
roi_result.formula = repmat("beta ~ -1 + cond + (1|SubjectID)", size(roi_result,1),1);
roi_result.AIC = repmat(AIC, size(roi_result,1),1);
roi_result.BIC = repmat(BIC, size(roi_result,1),1);
roi_result.group = repmat("OA", size(roi_result,1),1);
roi_result.model = repmat({roi_result_oa.model{1}}, size(roi_result,1),1);
condition_contrast_table = [condition_contrast_table; roi_result];

% PD
ContrastStats = GroupStats_pd.ttest(c);
roi_result = nirs.util.roiAverage(ContrastStats, ROI_ba_9_46, 'BA9_46');
AIC = roi_result_pd.model{1}.ModelCriterion.AIC;
BIC = roi_result_pd.model{1}.ModelCriterion.BIC;
roi_result.formula = repmat("beta ~ -1 + cond + (1|SubjectID)", size(roi_result,1),1);
roi_result.AIC = repmat(AIC, size(roi_result,1),1);
roi_result.BIC = repmat(BIC, size(roi_result,1),1);
roi_result.group = repmat("PD", size(roi_result,1),1);
roi_result.model = repmat({roi_result_pd.model{1}}, size(roi_result,1),1);
condition_contrast_table = [condition_contrast_table; roi_result];

disp('Main condition');
disp(main_condition_result);
disp('Contrast');
disp(condition_contrast_table);

%% PD HY stage contrast

% Remove cases with missing HY stage
selected_idx = ones(size(SubjStats_pd,2),1);
nan_idx = demographics_pd.hy == "NA";
selected_idx(nan_idx,1) = 0;
selected_idx = logical(selected_idx);
SubjStats_pd_hy = SubjStats_pd(selected_idx);
demographics_pd_hy = nirs.createDemographicsTable(SubjStats_pd_hy);
fprintf('Sum 1-2: %d\n', sum(demographics_pd_hy.hy == "HY_1_2"));
fprintf('Sum 3-4: %d\n', sum(demographics_pd_hy.hy == "HY_3_4"));

% Run group model
job = nirs.modules.MixedEffects();
job.formula = 'beta ~ -1 + cond:hy + (1|SubjectID)'; 
job.dummyCoding = 'full';
job.include_diagnostics = true;
GroupStats_pd_hy = job.run(SubjStats_pd_hy);

% Contrast stages for ST_walk
c = [0 0 -1 1 0 0];
ContrastStats = GroupStats_pd_hy.ttest(c);
roi_result_hy = nirs.util.roiAverage(ContrastStats, ROI_ba_9_46, 'BA9_46');
AIC = roi_result_pd.model{1}.ModelCriterion.AIC;
BIC = roi_result_pd.model{1}.ModelCriterion.BIC;
roi_result_hy.formula = repmat("beta ~ -1 + cond:hy + (1|SubjectID)", size(roi_result_hy,1),1);
roi_result_hy.AIC = repmat(AIC, size(roi_result_hy,1),1);
roi_result_hy.BIC = repmat(BIC, size(roi_result_hy,1),1);
roi_result_hy.group = repmat("PD", size(roi_result_hy,1),1);
roi_result_hy.model = repmat({roi_result_pd.model{1}}, size(roi_result_hy,1),1);

disp('HY stage contrast for ST_walk');
disp(roi_result_hy);

%% Interactions - YA

models_table = table();

job = nirs.modules.MixedEffects();
ya_formula = "beta ~ -1 + cond + hads_anxiety + cond:dt_cost_walk_speed + (1|SubjectID)"; 
job.formula = convertStringsToChars(ya_formula);
job.dummyCoding = 'full';
job.include_diagnostics = true;
GroupStats = job.run(SubjStats_ya);
roi_result_ya = nirs.util.roiAverage(GroupStats, ROI_ba_9_46, 'BA9_46');

disp('YA')
disp(roi_result_ya);
disp(roi_result_ya.model{1}.ModelCriterion.AIC);
disp(roi_result_ya.model{1}.ModelCriterion.BIC);

[AIC, BIC] = PlotDiagnostics(roi_result_ya.model{1}, "YA model", "Straight_walking_and_Aud_Stroop_dt_cost_walk_speed");

roi_result_ya.formula = repmat(ya_formula, size(roi_result_ya,1),1);
roi_result_ya.AIC = repmat(AIC, size(roi_result_ya,1),1);
roi_result_ya.BIC = repmat(AIC, size(roi_result_ya,1),1);
roi_result_ya.group = repmat("YA", size(roi_result_ya,1),1);

models_table = [models_table; roi_result_ya];

%% Interactions - OA

model_formulas = ["beta ~ -1 + cond + cond:balance + (1|SubjectID)",
    "beta ~ -1 + cond + cond:w12 + (1|SubjectID)",
    "beta ~ -1 + cond + cond:st_step_time_var + (1|SubjectID)",
    "beta ~ -1 + cond + cond:st_walk_speed + (1|SubjectID)",
    "beta ~ -1 + cond + cond:dt_cost_walk_speed + (1|SubjectID)"];

model_covariates = ["Straight_walking_balance", 
    "Straight_walking_w12", 
    "Straight_walking_st_step_time_var", 
    "Straight_walking_st_walk_speed", 
    "Straight_walking_and_Aud_Stroop_dt_cost_walk_speed"];

for idx=1:length(model_formulas)
    
    formula = convertStringsToChars(model_formulas(idx));
    covar = model_covariates(idx);

    job = nirs.modules.MixedEffects();
    job.dummyCoding = 'full';
    job.include_diagnostics = true;
    job.formula = formula;
    GroupStats = job.run(SubjStats_oa);
    roi_result_group = nirs.util.roiAverage(GroupStats, ROI_ba_9_46, 'BA9_46');

    [AIC, BIC] = PlotDiagnostics(roi_result_group.model{1}, "OA_model_" + covar , covar);
    
    roi_result_group.formula = repmat(model_formulas(idx), size(roi_result_group,1),1);
    roi_result_group.AIC = repmat(AIC, size(roi_result_group,1),1);
    roi_result_group.BIC = repmat(AIC, size(roi_result_group,1),1);
    roi_result_group.group = repmat("OA", size(roi_result_group,1),1);
    
    % Add to table
    models_table = [models_table; roi_result_group];
   
end

%% Interactions - PD

model_formulas = ["beta ~ -1 + cond + cond:balance + (1|SubjectID)",
    "beta ~ -1 + cond + cond:w12 + (1|SubjectID)",
    "beta ~ -1 + cond + cond:st_step_time_var + (1|SubjectID)",
    "beta ~ -1 + cond + cond:st_walk_speed + (1|SubjectID)",
    "beta ~ -1 + cond + cond:dt_cost_walk_speed + (1|SubjectID)",
    "beta ~ -1 + cond + cond:updrs_3_motor + (1|SubjectID)",];

model_covariates = ["Straight_walking_balance", 
    "Straight_walking_w12", 
    "Straight_walking_st_step_time_var", 
    "Straight_walking_st_walk_speed", 
    "Straight_walking_and_Aud_Stroop_dt_cost_walk_speed",
    "Straight_walking_updrs_3_motor"];

for idx=1:length(model_formulas)
    
    formula = convertStringsToChars(model_formulas(idx));
    covar = model_covariates(idx);
    
    job = nirs.modules.MixedEffects();
    job.dummyCoding = 'full';
    job.include_diagnostics = true;
    job.formula = formula;
    GroupStats = job.run(SubjStats_pd);
    roi_result_group = nirs.util.roiAverage(GroupStats, ROI_ba_9_46, 'BA9_46');

    [AIC, BIC] = PlotDiagnostics(roi_result_group.model{1}, "PD_model_" + covar, covar);
    
    roi_result_group.formula = repmat(model_formulas(idx), size(roi_result_group,1),1);
    roi_result_group.AIC = repmat(AIC, size(roi_result_group,1),1);
    roi_result_group.BIC = repmat(AIC, size(roi_result_group,1),1);
    roi_result_group.group = repmat("PD", size(roi_result_group,1),1);
    
    % Add to table
    models_table = [models_table; roi_result_group];
   
end

%% Pick variables of interest for hypotheses

hypothesis_table = table();

hypothesis_table = [hypothesis_table; ...
    main_condition_result(strcmp(main_condition_result.Contrast, "Straight_walking") & ~strcmp(main_condition_result.group, "YA"),:)];

hypothesis_table = [hypothesis_table; ...
    roi_result_hy(:,:)];

hypothesis_table = [hypothesis_table; ...
    condition_contrast_table(:,:)];

hypothesis_table = [hypothesis_table; ...
    models_table(strcmp(models_table.Contrast, "Straight_walking:balance"),:)];

hypothesis_table = [hypothesis_table; ...
    models_table(strcmp(models_table.Contrast, "Straight_walking:updrs_3_motor"),:)];

hypothesis_table = [hypothesis_table; ...
    models_table(strcmp(models_table.Contrast, "Straight_walking:w12"),:)];

hypothesis_table = [hypothesis_table; ...
    models_table(strcmp(models_table.Contrast, "Straight_walking:st_step_time_var"),:)];

hypothesis_table = [hypothesis_table; ...
    models_table(strcmp(models_table.Contrast, "Straight_walking:st_walk_speed"),:)];

hypothesis_table = [hypothesis_table; ...
    models_table(strcmp(models_table.Contrast, "Straight_walking_and_Aud_Stroop:dt_cost_walk_speed"),:)];

hypothesis_table = [hypothesis_table; ...
    models_table(strcmp(models_table.Contrast, "anxiety"),:)];

% Drop unneeded columns
hypothesis_table(:,[10,11]) = [];

disp(hypothesis_table);

%% Write table 

writetable(hypothesis_table, results_file)  

%% Draw

% Draw some figures.
GroupStats = GroupStats_ya;
GroupStats.probe.defaultdrawfcn='3D mesh (frontal)';
GroupStats.probe = GroupStats.probe.SetFiducialsVisibility(false);
GroupStats.draw('tstat', [-10 10], 'q < 0.05');

%% Print

% Print some figures.
GroupStats = GroupStats_oa;
GroupStats.probe.defaultdrawfcn='3D mesh (frontal)';
GroupStats.probe = GroupStats.probe.SetFiducialsVisibility(false);
folder = '../figures/group_activation_oa';
GroupStats.printAll('tstat', [-10 10], 'q < 0.05', folder, 'png')

%% Check assumptions

function [AIC, BIC] = PlotDiagnostics(model, titlestr, covar)

fitted_vals = model.Fitted;
residuals_raw = model.Residuals.Raw;
residuals_standardized = model.Residuals.Standardized;
leverage = model.Diagnostics.Leverage;

% residuals are (approximately) normally distributed (Q-Q plot)
figure;
subplot(3,2,1);
histfit(residuals_raw);
title('Histogram - raw residuals');
subplot(3,2,2);
qqplot(residuals_raw);
title('Q-Q - raw residuals');

% Residual vs fitted to assess possible heteroscedascity
% Plot fitted vs residuals
subplot(3,2,3);
h = plot(fitted_vals, residuals_standardized, 'bx');
xlabel("fitted");
ylabel("residual (standardized)");
hold on;

% Add 0 line and a polyline
p = polyfit(fitted_vals, residuals_standardized, 4);
px = linspace(min(fitted_vals), max(fitted_vals));
py = polyval(p, px);
plot(px, py, 'r-', 'LineWidth', 1);
xlim = [min(fitted_vals) max(fitted_vals)];
line(xlim,[0 0],'LineStyle','--', 'LineWidth', 1);
v = axis; % get current values
lo = min( v(1:2:end) ); % lower limit
up = max( v(2:2:end) ); % uppper limit
axis( [lo up lo up] );
hold off;
title('Fitted-Residual');

% Leverage vs residual
subplot(3,2,4);
plot(leverage, residuals_standardized, 'bx');
hold on;
xlabel("leverage");
ylabel("residual (standardized)");
p = polyfit(leverage, residuals_standardized, 4);
px = linspace(min(leverage), max(leverage));
py = polyval(p, px);
plot(px, py, 'r-', 'LineWidth', 1);
xlim = [min(leverage) max(leverage)];
line(xlim,[0 0],'LineStyle','--', 'LineWidth', 1);
hold off;
title('Leverage-Residual');

% Linearity in covariate
if covar ~= ""
    subplot(3,2,5);
    plot(residuals_raw, table2array(model.Variables(:,covar)), 'bx');
    title('Residual-Covariate');
end

sgtitle(titlestr) 

% Save the figure
fig = gcf;
set(gcf,'Position',[100 100 1000 700])
figname  = titlestr + "_" + covar + ".png";
figname = strrep(figname, " ", "_");
exportgraphics(fig,figname,'Resolution',300);

AIC = model.ModelCriterion.AIC;
BIC = model.ModelCriterion.BIC;

end