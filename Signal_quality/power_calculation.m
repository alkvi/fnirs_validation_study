% Get SCI power for additional movement artifact quantification. 
% Ran on R2021a, forked version
% https://github.com/alkvi/nirs-toolbox-fork/commits/stim_snirf_bids_fixes 
%% Load

% Load BIDS
my_data_dir = '../../Data/fNIRS_data/bids_dataset_snirf';
raw_data = nirs.bids.loadBIDS(my_data_dir, true);

% Also save the data to file
save('../../Data/mat_files/raw_data.mat','raw_data');

%% Preprocess 

% Label short-sep
job = nirs.modules.LabelShortSeperation();
job.max_distance = 8;
raw_data = job.run(raw_data);

% Trim baseline
job = nirs.modules.TrimBaseline();
job.preBaseline  = 1;
job.postBaseline = 1;
raw_data = job.run(raw_data);

% Set each duration to 20 seconds.
raw_data = nirs.design.change_stimulus_duration(raw_data,[],20);

job = nirs.modules.OpticalDensity();
od = job.run(raw_data);

%% Visualize

nirs.viz.TimeSeriesViewer(od);  

demographics = nirs.createDemographicsTable(od);
disp(demographics);

%% Set up common parameters

data_types = od(1).probe.link.type;
fcut = [0.5 2.5];
window = 5;
n_channels = 33;
fs = 10;

%% Power, setup 1

p1_idx = strcmp(demographics.session, "protocol1");
od_setup1 = od(p1_idx);

% Skip 118 which has only a few seconds of data
od_setup1(118) = [];

all_power_c1 = [];
all_power_c2 = [];
all_power_c3 = [];
for i=1:size(od_setup1,1)
    
    fprintf("On subject %d\n", i);
    
    od_data = od_setup1(i);
    time = od_data.time;
    input_data = od_data.data;
    stims_1 = od_data.stimulus('Straight_walking');
    stims_2 = od_data.stimulus('Stand_still_and_Aud_Stroop');
    stims_3 = od_data.stimulus('Straight_walking_and_Aud_Stroop');
    condition_data_1 = get_condition_data_channelwise(stims_1,time,input_data);
    condition_data_2 = get_condition_data_channelwise(stims_2,time,input_data);
    condition_data_3 = get_condition_data_channelwise(stims_3,time,input_data);

    data_in = condition_data_1;
    power_c1 = get_power_qt_nirs(data_in, data_types, fcut, window, n_channels, fs);
    power_c1 = reshape(power_c1,1,[]);
    data_in = condition_data_2;
    power_c2 = get_power_qt_nirs(data_in, data_types, fcut, window, n_channels, fs);
    power_c2 = reshape(power_c2,1,[]);
    data_in = condition_data_3;
    power_c3 = get_power_qt_nirs(data_in, data_types, fcut, window, n_channels, fs);
    power_c3 = reshape(power_c3,1,[]);
    
    all_power_c1 = [all_power_c1, power_c1];
    all_power_c2 = [all_power_c2, power_c2];
    all_power_c3 = [all_power_c3, power_c3];
    
    fprintf("Straight_walking - SCI Power - mean/std: %f/%f\n", mean(power_c1), std(power_c1));
    fprintf("Stand_still_and_Aud_Stroop - SCI Power - mean/std: %f/%f\n", mean(power_c2), std(power_c2));
    fprintf("Straight_walking_and_Aud_Stroop - SCI Power - mean/std: %f/%f\n", mean(power_c3), std(power_c3));

end

fprintf("\nTotal:\n");
fprintf("Straight_walking - SCI Power - mean/std: %f/%f\n", mean(all_power_c1), std(all_power_c1));
fprintf("Stand_still_and_Aud_Stroop - SCI Power - mean/std: %f/%f\n", mean(all_power_c2), std(all_power_c2));
fprintf("Straight_walking_and_Aud_Stroop - SCI Power - mean/std: %f/%f\n", mean(all_power_c3), std(all_power_c3));

%% Power, setup 2

p2_idx = strcmp(demographics.session, "protocol2");
od_setup2 = od(p2_idx);

all_power_c1 = [];
all_power_c2 = [];
for i=1:size(od_setup2,1)
    
    fprintf("On subject %d\n", i);
    
    od_data = od_setup2(i);
    time = od_data.time;
    input_data = od_data.data;
    stims_1 = od_data.stimulus('Straight_walking');
    stims_2 = od_data.stimulus('Navigated_walking');
    condition_data_1 = get_condition_data_channelwise(stims_1,time,input_data);
    condition_data_2 = get_condition_data_channelwise(stims_2,time,input_data);

    data_in = condition_data_1;
    power_c1 = get_power_qt_nirs(data_in, data_types, fcut, window, n_channels, fs);
    power_c1 = reshape(power_c1,1,[]);
    data_in = condition_data_2;
    power_c2 = get_power_qt_nirs(data_in, data_types, fcut, window, n_channels, fs);
    power_c2 = reshape(power_c2,1,[]);
    data_in = condition_data_3;
    
    all_power_c1 = [all_power_c1, power_c1];
    all_power_c2 = [all_power_c2, power_c2];
    
    fprintf("Straight_walking - SCI Power - mean/std: %f/%f\n", mean(power_c1), std(power_c1));
    fprintf("Navigated_walking - SCI Power - mean/std: %f/%f\n", mean(power_c2), std(power_c2));

end

fprintf("\nTotal:\n");
fprintf("Straight_walking - SCI Power - mean/std: %f/%f\n", mean(all_power_c1), std(all_power_c1));
fprintf("Navigated_walking - SCI Power - mean/std: %f/%f\n", mean(all_power_c2), std(all_power_c2));

%% Power, setup 3

p3_idx = strcmp(demographics.session, "protocol3");
od_setup3 = od(p3_idx);

all_power_c1 = [];
all_power_c2 = [];
for i=1:size(od_setup3,1)
    
    od_data = od_setup3(i);
    time = od_data.time;
    input_data = od_data.data;
    stims_1 = od_data.stimulus('Navigation');
    stims_2 = od_data.stimulus('Navigation_and_Aud_Stroop');
    condition_data_1 = get_condition_data_channelwise(stims_1,time,input_data);
    condition_data_2 = get_condition_data_channelwise(stims_2,time,input_data);

    data_in = condition_data_1;
    power_c1 = get_power_qt_nirs(data_in, data_types, fcut, window, n_channels, fs);
    power_c1 = reshape(power_c1,1,[]);
    data_in = condition_data_2;
    power_c2 = get_power_qt_nirs(data_in, data_types, fcut, window, n_channels, fs);
    power_c2 = reshape(power_c2,1,[]);
    data_in = condition_data_3;
    
    all_power_c1 = [all_power_c1, power_c1];
    all_power_c2 = [all_power_c2, power_c2];
    
    fprintf("On subject %d\n", i);
    fprintf("Navigation - SCI Power - mean/std: %f/%f\n", mean(power_c1), std(power_c1));
    fprintf("Navigation_and_Aud_Stroop - SCI Power - mean/std: %f/%f\n", mean(power_c2), std(power_c2));

end

fprintf("\nTotal:\n");
fprintf("Navigation - SCI Power - mean/std: %f/%f\n", mean(all_power_c1), std(all_power_c1));
fprintf("Navigation_and_Aud_Stroop - SCI Power - mean/std: %f/%f\n", mean(all_power_c2), std(all_power_c2));

%% Helper function

function condition_data = get_condition_data_channelwise(stims,time,data)

    condition_data = [];
    for onset_idx=1:length(stims.onset)
        onsets = stims.onset(onset_idx);
        durations = stims.dur(onset_idx);
        time_start_idx = find(time < onsets,1,'last');
        time_stop_idx = find(time > onsets+durations,1);
        interval_data = data(time_start_idx:time_stop_idx,:);
        condition_data = [condition_data; interval_data];
    end

end 
