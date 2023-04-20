

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% TempInt_EEG Phase Decoding on Perceptual Reports
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% In this script data are processed that resulted from an experiment in
% which participants reported one of two perceptual impressions of a visual
% stimulation after each trials. We use their concurrently recorded neural
% signals to predict the perceptual report in each trial from their neural
% activity before anything has been presented to them. Specifically, we
% extract the phase information from various neural frequencies, time
% points and electrodes for each participant and train a naive Bayes
% classifier with a subset of trials to then predict the perceptual reports
% of the remaining (non-trained) trials. To get an impression of how these
% classifier results would look like if our participants' perceptual
% reports were based on complete randomness, we additionally shuffle the
% assignment of trials to perceptual reports.

% Clean up the environment
clear all;
close all hidden;

% Define working directory
wd = 'D:\Experimental\TempInt_EEG\';

% Load the Cosmo functions
addpath([wd, 'functions\CoSMoMVPA-master\']);
addpath(genpath([wd, 'functions\CoSMoMVPA-master\']));

% Load the fieldtrip functions
addpath([wd, 'functions\fieldtrip-20190419\']);
ft_defaults;


%% Define general variables

% Define subject numbers
vpn = [2,3,6,7,12,14,15,16,19,22,26,27,28,30,33];

% Define the channel variable
chans = 1:64;

% Define the permutation number
perm_n = 1000;

% Load channel information
load([wd, 'results\chanlocs.mat']);

% Load participant data to prepare general information
load([wd, 'results\tfa\fourierspec\ref_mast\v2_fourierspec.mat'], ...
    'tfdata', 'times', 'freqs', 'trlinfo');

% Predefine the resulting matrix
cr_perm = zeros(perm_n, 720);


%% Prepare the classifier structure

% Prepare general information that are the same for all participants
ds_tl.fa.chan = ...
    reshape(permute(repmat((1:length(chans))', ...
    [1,length(freqs),length(times)]), [1,2,3]), ...
    [1,length(chans)*length(freqs)*length(times)]);
ds_tl.fa.freq = ...
    reshape(permute(repmat((1:length(freqs))', ...
    [1,length(chans),length(times)]), [2,1,3]), ...
    [1,length(chans)*length(freqs)*length(times)]);
ds_tl.fa.time = ...
    reshape(permute(repmat((1:length(times))', ...
    [1,length(chans),length(freqs)]), [2,3,1]), ...
    [1,length(chans)*length(freqs)*length(times)]);

% Prepare parameter labels on which the classifier will be trained
ds_tl.a.fdim.labels = {'chan', 'freq', 'time'};
ds_tl.a.fdim.values = {{chanloc.labels}, freqs, times};
ds_tl.a.meeg.samples_field = 'phasespctrm';

% Define centers on which the searchlight should act
time_ind = 3:5:48;
freq_ind = 6:5:31;
chan_ind = [3,10,18,26,30,32,36,37,45,47,55,63];

% Define indices for those channel-time-frequency combination that will be
% considered in the classifier
center_ids_mat = find(ismember(ds_tl.fa.chan, chan_ind) & ...
    ismember(ds_tl.fa.freq, freq_ind) & ismember(ds_tl.fa.time, time_ind));


%% Fill example data into the classifier structure

% Berechnung der Phasen-Information
phang = angle(tfdata);

% Neu-Strukturierung der Phasen-Daten
phang_reshaped = reshape(permute(phang, [4,1,3,2]), ...
    [size(phang,4), size(phang,1)*size(phang,2)*size(phang,3)]);

% Save data in structure that will be forwarded to the classifier
ds_tl.samples = phang_reshaped;

% Define the neighborhood for each dimensions (neighbors are not
% incorporated)
chan_nbrhood = cosmo_meeg_chan_neighborhood(ds_tl, ...
    'count', 9, 'chantype', 'eeg');
freq_nbrhood = cosmo_interval_neighborhood(ds_tl, 'freq', ...
    'radius', 2);
time_nbrhood = cosmo_interval_neighborhood(ds_tl, 'time', ...
    'radius', 2);

% Cross neighborhoods for chan-time-freq searchlight
nbrhood = cosmo_cross_neighborhood(ds_tl, ...
    {chan_nbrhood, freq_nbrhood, time_nbrhood});


%% Use classifier on each subject's data

% Participant loop
for v = 1:length(vpn)

    
    %% Prepare participant data
    
    % Load participant data to prepare general information
    load([wd, 'results\tfa\fourierspec\ref_mast\v', num2str(vpn(v)), ...
        '_fourierspec.mat']);

    % Calculate phase information
    phang = angle(tfdata);

    % Restructure phase data
    phang_reshaped = reshape(permute(phang, [4,1,3,2]), ...
        [size(phang,4), size(phang,1)*size(phang,2)*size(phang,3)]);

    % Permutation loop
    for n = 1:perm_n

        % Randomize randomness
        rng('shuffle');
        
        % Randomize reports to each trial
        trlinfo = trlinfo(randperm(length(trlinfo)));
        
        % Use the randomized trial information to add phase information to
        % the classifier structure
        trlinfo(trlinfo == 3) = 2;
        [trlinfo_sorted, sort_ind] = sort(trlinfo, 'ascend');
        ds_tl.samples = phang_reshaped(sort_ind,:);
        trl_ind = [1:length(sort_ind)/2, 1:length(sort_ind)/2]';

        % Add trlinfo to the structure
        ds_tl.sa.trialinfo = [trlinfo_sorted', trl_ind];

        % Add a matrix with trials*freqs
        ds_tl.sa.cumtapent = ones(length(trl_ind),1);

        % Set targets (trial condition)
        ds_tl.sa.targets = ds_tl.sa.trialinfo(:,1);

        % Set initial chunks
        ds_tl.sa.chunks = ds_tl.sa.trialinfo(:,2);


        %% Set up classifier parameters

        % Reset chunks: use four chunks
        nchunks = 4;
        ds_tl.sa.chunks = cosmo_chunkize(ds_tl, nchunks);

        % Prepare a take-one-fold out cross validation
        partitions = cosmo_nchoosek_partitioner(ds_tl, 1);
        partitions = cosmo_balance_partitions(partitions, ds_tl);

        % Define the measure and its arguments
        measure = @cosmo_crossvalidation_measure;
        measure_args = struct();
        measure_args.classifier = @cosmo_classify_naive_bayes;
        measure_args.partitions = partitions;

        % Combine measurement arguments to one structure
        measure_args = cosmo_structjoin(measure_args, ...
            {'average_train_count', 5, ...
            'average_train_resamplings', 3});


        %% Run the classifier

        % Classifier is only executed for center_ids
        sl_map = cosmo_searchlight(ds_tl, nbrhood, ...
            measure, measure_args, 'center_ids', center_ids_mat);

        % Save the classifier result
        cr_perm(n,:) = sl_map.samples;
        
    end
    % Permutation loop
    
    % Save classification results of the permuted datasets
    save([wd, 'results\behavioral\phase_reverse\mvpa\mvpa_bayes_', ...
        'searchlight_whole_perm_v', num2str(vpn(v)), '.mat'], ...
        'cr_perm', 'times', 'freqs', 'chans', 'time_ind', 'freq_ind', ...
        'chan_ind');
    
end
% Participant loop


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% TempInt_EEG Phase Decoding on Perceptual Reports
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

