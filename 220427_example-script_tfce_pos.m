

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Start: Conduct inference statistics according to TFCE
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Clear both workspace and open figure windows
clear all;
close all hidden;


% Define the working directory
wd = '/home/matlab/Alex/TempInt_InfStatComp/';


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Load functions
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Load parallel package
pkg load parallel

% Load image package
pkg load image

% Consider additional functions
addpath([wd, 'functions/']);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Define comparison parameters
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Define the study
studlab = 'TempInt_EEG';

% Define conditions that shall be compared
complab = 'ZiMvsZvM';

% Define condition inidacators
condind = [1,3];


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Load time-frequency data
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Define participant variable
vpn = [2,3,6,7,12,14,15,16,19,22,26,27,28,30,33];

% Participant loop
for v = 1:length(vpn)
    
    % Load data
    load([wd, 'results/', studlab, '/tfa/wavelet/', ...
        'v', num2str(vpn(v)), '_fourierspec.mat'], ...
        'tfdata', 'chans', 'times', 'freqs', 'trlinfo');
    
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
    %% Extract trials for the condition comparison
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
    
    
    % Calculate the phase angle
    phase_mat = exp(1i .* angle(tfdata));
    
    % Define an indicator for trials of interest
    trl_ind = trlinfo == condind(1) | trlinfo == condind(2);
    
    % Reduce the dataset to trials of interest and time points after
    % target onset
    data_cell{v} = phase_mat(:,times < 0,:,trl_ind);
    
    % Reduce the trial information to trials of interest and drop the
    % running number
    trlinfo_cell{v} = trlinfo(trl_ind);
    
    % Clean up
    clear('tfdata', 'phase_mat', 'trlinfo');
    
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
    %% Display loading progression
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
    
    
    % Display the participant number that was successfully loaded
    disp(['Participant ', num2str(v), '/', num2str(length(vpn)), ...
        ' was loaded.']);
    
end
% Participant Loop

% Also restrict the time variable to time points before target onset
times = times(times < 0);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Define a topographical channel structure
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Load a variable holding information on the channel neighbourhood
nghbrhd = getfield(load([wd, 'results/', ...
    'clustering_fundamentals.mat'], ...
    'chantopo_num'),'chantopo_num');

% Hold information on channels that occur multiple times in the
% neighbourhood
dropchan(1,:) = getfield(load([wd, 'results/', ...
    'clustering_fundamentals.mat'], ...
    'dropchan_x'),'dropchan_x');
dropchan(2,:) = getfield(load([wd, 'results/', ...
    'clustering_fundamentals.mat'], ...
    'dropchan_y'),'dropchan_y');


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Define permutation parameters
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Anzahl zu berechnender Permutationen
rp_n = 1000;


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Define clustering parameters
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Define the theoretical connectivity matrix
conn = zeros(3,3,3,3);
conn(2,2,2,2) = 1;
conn(1,2,2,2) = 1; conn(3,2,2,2) = 1;
conn(2,1,2,2) = 1; conn(2,3,2,2) = 1;
conn(2,2,1,2) = 1; conn(2,2,3,2) = 1;
conn(2,2,2,1) = 1; conn(2,2,2,3) = 1;


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Define parallel processing parameters
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Define the number of cores
cores_n = 40;

% Define the number of permutations each core needs to process
perm_n = rp_n / cores_n;

% Define a core indicator so that each core processes different information
core_ind = 1:cores_n;


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Calculate the true dataset
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Try to load the true data, otherwise calculate them anew
try

    % Display progress
    disp(['Trying to load true data.']);

    % Load true data
    load([wd, 'results/', studlab, '/tfa/pos/surrogate/pos_tm.mat'], ...
        'data_tm');
        
    % Display progress
    disp(['True data have been calculated before and are now loaded.']);
    
catch
    
    % Display progress
    disp(['True data have not been calculated before. ', ...
        'They are now calculated.']);
    
    % Participant loop
    for v = 1:length(vpn)
        
        % Calculate the itc
        itc_tm(v,1,:,:,:) = ...
            abs(mean(data_cell{v}(:,:,:,trlinfo_cell{v} == ...
                condind(1)),4));
        itc_tm(v,2,:,:,:) = ...
            abs(mean(data_cell{v}(:,:,:,trlinfo_cell{v} == ...
                condind(2)),4));
        itc_tm(v,3,:,:,:) = ...
            abs(mean(data_cell{v},4));
        
        % Calculate the pos
        data_tm(v,:,:,:) = squeeze(itc_tm(v,1,:,:,:) + ...
                                   itc_tm(v,2,:,:,:) - ...
                              2 .* itc_tm(v,3,:,:,:));
                              
        % Display participant number that was successfully processed
        disp(['Participant ', num2str(v), '/', num2str(length(vpn)), ...
            ' was processed.']);
        
    end
    % Participant loop
    
    % Display progress
    disp(['True data have been calculated are now being saved.']);
    
    % Save the true data
    save('-mat7-binary', [wd, 'results/', studlab, '/tfa/pos/', ...
        'surrogate/pos_tm.mat'], ...
        'data_tm', 'itc_tm', 'chans', 'times', 'freqs');
        
    % Display progress
    disp(['True data have been saved.']);

end
% Try to load the true data, otherwise calculate them anew


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Generate a surrogate dataset
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Define the permutations as a function to run them in parallel on 
% multiple cores
function  data_pmp = permfun(u, perm_n, vpn, chans, times, freqs, ...
    condind, trlinfo_cell, data_cell)
    
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
    %% Predefine resulting matrices
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
    
    
    % Predefine the matrix in which each permutation result will be saved
    data_pmp = zeros(perm_n, length(vpn), length(chans), length(times), ...
        length(freqs));
    
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
    %% Start permutations
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
    
    
    % Start permutation loop
    for p = 1:perm_n
        
        
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
        %% Data Permutation
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
        
        
        % Participant loop
        for v = 1:length(vpn)
            
            % Randomize randomness
            rand('state', 'reset');
            
            
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
            %% Randomize trials
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
            
                            
            % Define a vector of random indicators
            randind = randperm(length(trlinfo_cell{v}));

            % Save the randomized order of trials in the trial information
            trlinfo_cell{v} = trlinfo_cell{v}(randind);
            
            
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
            %% Calculate POS
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
            
            
            % Calculate the itc
            itc_pm(1,:,:,:) = ...
                abs(mean(data_cell{v}(:,:,:,trlinfo_cell{v} == ...
                    condind(1)),4));
            itc_pm(2,:,:,:) = ...
                abs(mean(data_cell{v}(:,:,:,trlinfo_cell{v} == ...
                    condind(2)),4));
            itc_pm(3,:,:,:) = ...
                abs(mean(data_cell{v},4));
            
            % Calculate the pos
            data_pmp(p,v,:,:,:) = squeeze(itc_pm(1,:,:,:) + ...
                                          itc_pm(2,:,:,:) - ...
                                     2 .* itc_pm(3,:,:,:));
            
        end
        % Versuchspersonen-Schleife
        
        
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
        %% Display permutation progression
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


        % Display progress
        if u == 1
            disp(['Permutation ', num2str(p), ' has been processed.']);
        end

    end
    % Permutation loop

endfunction
% Define the permutations as a function to run them in parallel on
% multiple cores

% Add a variable to the function that is defined for each core
parapermfun = @(u)permfun(u, perm_n, vpn, chans, times, freqs, ...
    condind, trlinfo_cell, data_cell);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Start the permutation procedure
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Only perform the permutations if there is no surrogate dataset to be
% loaded
try

    % Display loading progression
    disp(['Permutations are trying to be loaded.']);

    % Predefine the resulting matrix
    data_pm = zeros(rp_n, length(vpn), length(chans), length(times), ...
        length(freqs));

    % Participant loop
    for v = 1:length(vpn)
        
        % Load data separately for each participant
        load([wd, 'results/', studlab, '/tfa/pos/surrogate/', ...
            'pos_pm_v', num2str(vpn(v)), '.mat'], ...
            'data_pmp');
            
        % Save all participant data into the same matrix
        data_pm(:,v,:,:,:) = data_pmp;
        
        % Display loading progression
        disp(['Participant ', num2str(v) ' was loaded.']);
            
    end
    % Participant loop

catch

    % Display loading progression
    disp(['No permutation data could be loaded.', ...
        ' Permutations start now.']);

    % Start timer
    tic

    % Start permutations
    data_pmp_cell = pararrayfun(cores_n, parapermfun, core_ind, ...
        'UniformOutput', false);

    % End timer
    toc
    
    % Display processing progression
    disp(['Permutations have been finished! Cells are now translated', ...
        ' into a matrix.']);


    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
    %% Translate cell into matrix
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


    % Core Loop
    for u = 1:length(data_pmp_cell)

        % Core Statement
        if u == 1

            % Create a matrix for the permuted data
            data_pm = data_pmp_cell{u};

        else

            % Add further permuted data to the matrix
            data_pm = cat(1, data_pm, data_pmp_cell{u});

        end
        % Core Statement

    end
    % Core Loop
    
    % Calculate the mean of all permutations
    data_pm_avg = squeeze(mean(data_pm,1));

    % Calculate the standard deviation of all permutations
    data_pm_std = squeeze(std(data_pm,[],1));
    
    % Display processing progression
    disp(['Translations have been finished! Data is saved now.']);


    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
    %% Save permuted data
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

        
    % Save the mean and standard deviation of all permutations
    save('-mat7-binary', [wd, 'results/', studlab, '/tfa/pos/surrogate/', ...
        '/pos_pm_avgstd.mat'], ...
        'data_pm_avg', 'data_pm_std');
        
        
    % Participant loop
    for v = 1:length(vpn)
        
        % Extract trial data
        data_pmp = squeeze(data_pm(:,v,:,:,:));
        
        % Save data separately for each participant
        save('-mat7-binary', [wd, 'results/', studlab, '/tfa/pos/', ...
            'surrogate/pos_pm_v', num2str(vpn(v)), '.mat'], ...
            'data_pmp');
            
        % Display saving progression
        disp(['Participant ', num2str(v), '/', num2str(length(vpn)), ...
            ' has been saved.']);
            
    end
    % Participant loop

end
% Only perform the permutations if there is no surrogate dataset to be
% loaded


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Calculate z statistics
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Perform a t-test against zero
data_tma = squeeze(mean(data_tm,1));

% Perform a t-test against zero
data_pma = squeeze(mean(data_pm,2));

% Additionally, determine the average for the permutation dimension
data_pma_avg = squeeze(mean(data_pma,1));

% Additionally, determine the standard deviation for the permutation
% dimension
data_pma_std = squeeze(std(data_pma,[],1));

% z transform true data using the mean and standard deviation of the
% permutation distribution and calculate the participant mean
data_tmaz = (data_tma - data_pma_avg) ./ data_pma_std;

% z transform permuted data using the mean and standard deviation of the
% permutation distribution
data_pmaz = (data_pma - ...
    permute(repmat(data_pma_avg, [1,1,1,rp_n]),[4,1,2,3])) ./ ...
    permute(repmat(data_pma_std, [1,1,1,rp_n]),[4,1,2,3]);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Organize channel information topographically
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Topographical x loop
for x = 1:size(nghbrhd,1)

    % Topographical y loop
    for y = 1:size(nghbrhd,2)

        % Check whether the topographical spot is a channel
        if nghbrhd(x,y) > 0

            % Copy data at the right topographical spot of the
            % true data
            topo_tmaz(x,y,:,:) = ...
                squeeze(data_tmaz(nghbrhd(x,y),:,:));
            
            % Copy data at the right topographical spot of the
            % permuted data
            topo_pmaz(:,x,y,:,:) = ...
                squeeze(data_pmaz(:,nghbrhd(x,y),:,:));

        end
        % Check whether the topographical spot is a channel

    end
    % Topographical Y Loop

end
% Topographical X Loop

% Save the z transformed topographical data
save('-mat7-binary', [wd, 'results/', studlab, '/tfa/pos/', ...
    'surrogate/pos_topo_zval.mat'], ...
    'topo_tmaz', 'data_tmaz', 'data_tma');


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Apply the threshold-free cluster-enhancement to the true data
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Define the integral steps
steps = .1:.1:5.0;

% Predefine the final matrix that contains the tfce values
tfce_tmt = zeros(size(nghbrhd,1), size(nghbrhd,2), length(times), ...
    length(freqs));

% Integral Step Loop
for h = 1:length(steps)

    % Keep all significant positive data points
    topo_sig = topo_tmaz > steps(h);

    % Cluster all potential cluster elements
    topo_conn = bwlabeln(topo_sig, conn);
    
    % Correct for channels that occured double
    for i = 1:size(dropchan,2)
        topo_conn(dropchan(1,i),dropchan(2,i),:) = 0;
    end

    % Extract the number of cluster elements in each cluster
    topo_reg = regionprops(topo_conn, 'Area');

    % Data Point Loop
    for i = 1:numel(topo_conn)

        % Perform the following steps only if the data point is still
        % greater than the intergral step
        if topo_conn(i) > 0

            % Add the tfce results for this step size to each data point
            tfce_tmt(i) = tfce_tmt(i) + ...
                sqrt(topo_reg(topo_conn(i)).Area) * steps(h)^2;

        end
        % Perform the previous steps only if the data point is still
        % greater then the intergral step

    end
    % Data Point Loop

end
% Integral Step Loop


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Reshape the matrix back to a normal channel structure
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Duplicate the neighbourhood matrix
nghbrhd_red = nghbrhd;

% Reduce the neighbourhood matrix of double channel occurences
for i = 1:size(dropchan,2)
    nghbrhd_red(dropchan(1,i), dropchan(2,i)) = 0;
end

% Predefine the resulting matrix
tfce_tm = zeros(length(chans), length(times), length(freqs));

% Topographical x loop
for x = 1:size(nghbrhd_red,1)

    % Topographical y loop
    for y = 1:size(nghbrhd_red,2)
    
        % Only copy data if the spatial position corresponds to a channel
        if nghbrhd_red(x,y) > 0
        
            % Copy data from a topographical to a normal channel
            % structure
            tfce_tm(nghbrhd_red(x,y),:,:) = squeeze(tfce_tmt(x,y,:,:));
            
        end
        % Only copy data if the spatial position corresponds to a channel
    
    end
    % Topographical y loop
    
end
% Topographical x loop


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Save the cluster results of the true data
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Save the cluster results of the true data
save('-mat7-binary', [wd, 'results/', studlab, '/tfa/pos/global/', ...
    'tfce/pos_tm_clust.mat'], ...
    'tfce_tm', 'tfce_tmt');


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Apply the threshold-free cluster-enhancement to the permuted data
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Define the tfce of the permuted data as a function to run them in
% parallel on multiple cores
function [tfce_pm_temp] = permfun(u, perm_n, times, chans, freqs, ...
    topo_pmaz, conn, nghbrhd, dropchan)
    
    
    % Define the integral steps
    steps = .1:.1:5.0;
    
    % Define the final matrix containing the largest tfce results of all
    % permutations
    tfce_pm_temp = zeros(1, perm_n);
    
    % Define the permutations that need to be processed by this core
    perm_ind = (1:perm_n) + ((u-1) * perm_n);
    
    % Permutation loop
    for p = 1:perm_n

        % Extract data of one permutation
        topo_pmazp = squeeze(topo_pmaz(perm_ind(p),:,:,:,:));
        
        % Create a matrix containing the tfce results for a permutation
        tfce_pm_mat = zeros(size(nghbrhd,1), size(nghbrhd,2), ...
            length(times), length(freqs));

        % Integral Step Loop
        for h = 1:length(steps)

            % Keep all significant positive data points
            topo_sig = topo_pmazp > steps(h);

            % Cluster all potential cluster elements
            topo_conn = bwlabeln(topo_sig, conn);
            
            % Correct for channels that occured double
            for i = 1:size(dropchan,2)
                topo_conn(dropchan(1,i),dropchan(2,i),:) = 0;
            end

            % Extract the number of cluster elements in each
            % cluster
            topo_reg = regionprops(topo_conn, 'Area');

            % Data Point Loop
            for i = 1:numel(topo_conn)

                % Perform the following steps only if the data point
                % is still greater than the intergral step
                if topo_conn(i) > 0

                    % Add the tfce results for this step size to each
                    % data point
                    tfce_pm_mat(i) = tfce_pm_mat(i) + ...
                        sqrt(topo_reg(topo_conn(i)).Area) * ...
                        steps(h)^2;

                end
                % Perform the previous steps only if the data point
                % is still greater then the intergral step

            end
            % Data Point Loop

        end
        % Integral Step Loop

        % Extract the largest tfce result from this channel
        tfce_pm_temp(1,p) = max(tfce_pm_mat(:));
        
        % Display progress
        if u == 1
            disp(['Permutation ', num2str(p), '/', num2str(perm_n), ...
                ' has been processed.']);
        end
        % Display progress        

    end
    % Permutation Loop
    
endfunction
% Define the permutations as a function to run them in parallel on
% multiple cores

% Add a variable to the function that is defined for each core
parapermfun = @(u)permfun(u, perm_n, times, chans, freqs, topo_pmaz, ...
    conn, nghbrhd, dropchan);

% Start permutations
[tfce_pm_cell] = ...
    pararrayfun(cores_n, parapermfun, core_ind, 'UniformOutput', false);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Translate cells into matrix
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Core Loop
for u = 1:length(tfce_pm_cell)

    % Core Statement
    if u == 1

        % Create a matrix for the permuted data
        tfce_pm = tfce_pm_cell{u};
        
    else

        % Add further permuted data to the matrix
        tfce_pm = [tfce_pm, tfce_pm_cell{u}];
        
    end
    % Core Statement

end
% Core Loop


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Save the cluster results of the permuted data
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Save the cluster results of the permuted data
save('-mat7-binary', [wd, 'results/', studlab, '/tfa/pos/global/', ...
    'tfce/pos_pm_clust.mat'], ...
    'tfce_pm');


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Start: Conduct inference statistics according to TFCE
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

