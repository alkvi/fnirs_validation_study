% Testing
% Ran on R2021a, forked version of NIRS toolbox
% https://github.com/alkvi/nirs-toolbox-fork/tree/phd_study_2
%% Load

% Load the NIRx probe used
probe_folder = "../Data/fNIRS_data/nirx_probe";
nirx_probe = nirs.io.loadNIRxProbe(probe_folder,true);

% Load BIDS
my_data_dir = '../Data/fNIRS_data/bids_dataset_snirf';
raw_data = nirs.bids.loadBIDS(my_data_dir, true, nirx_probe);

% Also save the data to file
save('../Data/mat_files/raw_data.mat','raw_data');

% ..or load
%load('../Data/mat_files/raw_data.mat');

demographics = nirs.createDemographicsTable(raw_data);
disp(demographics);

%% Add age to demographics

age_data = readtable('../Data/basic_demographics.csv'); 
demographics.age = NaN(height(demographics),1);

% Add age data
for idx=1:height(age_data)
    subj_id_seek = string(age_data.subject(idx)); 
    match_idx = strcmp(demographics.SubjectID, subj_id_seek);
    if sum(match_idx) < 1
        continue
    end
    demographics(match_idx,:).age = ones(sum(match_idx),1) * age_data(idx,:).age;
end

job = nirs.modules.AddDemographics;
job.demoTable = demographics;
job.varToMatch = 'UUID';
raw_data = job.run(raw_data);

%% Prune bad channels in data

bad_ch_file = "../Data/temp_data/sci_bad_channels.csv";
bad_ch_table = readtable(bad_ch_file);

% Loop through bad channel table for each src-det pair.
for i = 1:size(bad_ch_table,1)
    link = bad_ch_table(i,:);
    
    % Find which subjects have this pair marked as bad.
    subject_list = string(bad_ch_table(i,:).subject).split(',');
    if strcmp(subject_list,"None")
        continue
    end
    
    % raw_data(j).data is TxN where N is number of channels
    % replace raw_data(j).data with NaN in col N for bad channel in N
    for j = 1:length(subject_list)
        subj_parts = subject_list(j).split('-');
        subj = subj_parts(1);
        protocol = subj_parts(2);
        subj_idx = find(strcmp(demographics.SubjectID,subj) & strcmp(demographics.session,protocol));
        idx_to_null = raw_data(subj_idx).probe.link.source == link.src ...
            & raw_data(subj_idx).probe.link.detector == link.det;
        
        fprintf("Nulling src-det %d-%d for subj %s\n", link.src, link.det, subj);
        raw_data(subj_idx).data(:,idx_to_null) = NaN;
    end
end

%% Preprocess 

% Label short-sep
job = nirs.modules.LabelShortSeperation();
job.max_distance = 8;
raw_data = job.run(raw_data);

% Extra short-sep labelling:
% Some short-detectors (virtual detector 8-15) make up links 
% with not just the source they're attached to, but ones further away.
% these are still short-separation (so e.g. 1-14 and 1-15 detect the same
% signal). Mark these also. NOTE: the column is called "ShortSeperation" (sic).
for i = 1:numel(raw_data)
    short_separation_idx =(raw_data(i).probe.link.detector >= 8);
    raw_data(i).probe.link.ShortSeperation(short_separation_idx) = 1;
end

% raw_data(4) had an extra onset.
raw_data(4).stimulus('Straight_walking').onset(1) = [];

% Trim baseline
job = nirs.modules.TrimBaseline();
job.preBaseline  = 1;
job.postBaseline = 1;
raw_data = job.run(raw_data);

% Set each duration to 20 seconds.
raw_data = nirs.design.change_stimulus_duration(raw_data,[],20);

job = nirs.modules.OpticalDensity();
od = job.run(raw_data);

job = nirs.modules.BeerLambertLaw();
hb = job.run(od);

%% Visualize

nirs.viz.TimeSeriesViewer(hb);  

demographics = nirs.createDemographicsTable(hb);
disp(demographics);

figure;
raw_data(1).probe.defaultdrawfcn='3D mesh';
raw_data(1).probe.draw;

%% Get combined measures

% Add CBSI
job = nirs.modules.CalculateCBSI();
hb = job.run(hb);

% Add HbT
job = nirs.modules.CalculateTotalHb();
hb = job.run(hb);

% Add HbDiff
job = nirs.modules.CalculateHbDiff();
hb = job.run(hb);

job=nirs.modules.KeepTypes;
job.types={'cbsi', 'hbo', 'hbr', 'hbt', 'hbd'};
hb=job.run(hb);

%% Run on protocol 1

p1_idx = strcmp(demographics.session, "protocol1");
hb_setup1 = hb(p1_idx);

% Skip 118 which has only a few seconds of data
% Exclude also 84, 104 who did not stand still during rest.
hb_setup1([118,84,104]) = [];

nirs.viz.TimeSeriesViewer(hb_setup1);  

%% Then run GLM 

% Only CBSI
hb_setup1_cbsi = hb_setup1;
job=nirs.modules.KeepTypes;
job.types={'cbsi'};
hb_setup1_cbsi=job.run(hb_setup1_cbsi);
job = nirs.modules.GLM();
job.type='AR-IRLS';
job.AddShortSepRegressors = true;
SubjStats_cbsi = job.run(hb_setup1_cbsi); 
save('../Data/mat_files/SubjStats_setup_1_cbsi_vardpf.mat','SubjStats_cbsi');
clear SubjStats_cbsi;
clear hb_setup1_cbsi;

% Only HbO
hb_setup1_hbo = hb_setup1;
job=nirs.modules.KeepTypes;
job.types={'hbo'};
hb_setup1_hbo=job.run(hb_setup1_hbo);
job = nirs.modules.GLM();
job.type='AR-IRLS';
job.AddShortSepRegressors = true;
SubjStats_hbo = job.run(hb_setup1_hbo); 
save('../Data/mat_files/SubjStats_setup_1_hbo_vardpf.mat','SubjStats_hbo');
clear SubjStats_hbo;
clear hb_setup1_hbo;

% Only HHb
hb_setup1_hbr = hb_setup1;
job=nirs.modules.KeepTypes;
job.types={'hbr'};
hb_setup1_hbr=job.run(hb_setup1_hbr);
job = nirs.modules.GLM();
job.type='AR-IRLS';
job.AddShortSepRegressors = true;
SubjStats_hbr = job.run(hb_setup1_hbr); 
save('../Data/mat_files/SubjStats_setup_1_hbr_vardpf.mat','SubjStats_hbr');
clear SubjStats_hbr;
clear hb_setup1_hbr;

% Only HbT
hb_setup1_hbt = hb_setup1;
job=nirs.modules.KeepTypes;
job.types={'hbt'};
hb_setup1_hbt=job.run(hb_setup1_hbt);
job = nirs.modules.GLM();
job.type='AR-IRLS';
job.AddShortSepRegressors = true;
SubjStats_hbt = job.run(hb_setup1_hbt); 
save('../Data/mat_files/SubjStats_setup_1_hbt_vardpf.mat','SubjStats_hbt');
clear SubjStats_hbt;
clear hb_setup1_hbt;

% Only HbD
hb_setup1_hbd = hb_setup1;
job=nirs.modules.KeepTypes;
job.types={'hbd'};
hb_setup1_hbd=job.run(hb_setup1_hbd);
job = nirs.modules.GLM();
job.type='AR-IRLS';
job.AddShortSepRegressors = true;
SubjStats_hbd = job.run(hb_setup1_hbd); 
save('../Data/mat_files/SubjStats_setup_1_hbd_vardpf.mat','SubjStats_hbd');
clear SubjStats_hbd;
clear hb_setup1_hbd;

