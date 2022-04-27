

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Start: Conduct inference statistics according to Maris & Oostenveld
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
%% Load ERP data
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Define participant variable
vpn = [2,3,6,7,12,14,15,16,19,22,26,27,28,30,33];

% Participant loop
for v = 1:length(vpn)
    
    % Load data
    load([wd, 'results/', studlab, '/tfa/wavelet/', ...
        'v', num2str(vpn(v)), '_fourierspec.mat'], ...
        'tfdata', 'chans', 'times', 'freqs', 'trlinfo');
    
    % Calculate the power
    power_mat = tfdata .* conj(tfdata);
    
    % Clear space
    clear('tfdata');
    
    % Calculate the trial average for each condition and participant
    data_tm(v,1,:,:,:) = mean(power_mat(:,:,:,trlinfo == condind(1)),4);
    data_tm(v,2,:,:,:) = mean(power_mat(:,:,:,trlinfo == condind(2)),4);
    
    % Clear space
    clear('power_mat');    
    
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
    %% Display loading progression
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
    
    
    % Display the participant number that was successfully loaded
    disp(['Participant ', num2str(v), '/', num2str(length(vpn)), ...
        ' was loaded.']);
    
end
% Participant Loop


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
cores_n = 32;

% Define a core indicator so that each core processes different information
cores_ind = 1:cores_n;


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Baseline normalize power data
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Calculate the baseline
bl = repmat(mean(data_tm,4), [1,1,1,size(data_tm,4),1]);
    
% Baseline the true data
power_tm = 100 .* ((data_tm - bl) ./ bl);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Save the true dataset
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Save the true data
save('-mat7-binary', [wd, 'results/', studlab, '/tfa/power/', ...
    'surrogate/spec_bl/power_tm.mat'], ...
    'power_tm', 'chans', 'times');


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Prepare the full permutation procedure
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Define a matrix that systematically swaps averages between true and
% permuted data
combmat = permute(repmat(allcomb( ...
    [1,2],[1,2],[1,2],[1,2],[1,2],[1,2],[1,2],[1,2],[1,2],[1,2], ...
    [1,2],[1,2],[1,2],[1,2],[1,2]),[1,1,2]),[2,3,1]);
combmat(:,2,:) = abs(combmat(:,2,:)-3);

% As soon as there is a 2 in the first column, averages need to be swapped
shiftmat = squeeze(combmat(:,1,:)) == 2;
    
% Clear things up
clear('combmat');

% Dissolve the channel, time and frequency dimension into one
data_perm = reshape(data_tm, [size(data_tm,1), size(data_tm,2), ...
    size(data_tm,3) * size(data_tm,4) * size(data_tm,5)]);
    

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Create permutation indicators for each core
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Define the number of permutations each core needs to process
core_perm_n = size(shiftmat,2) / cores_n;

% Core Loop
for c = 1:cores_n
    
    % Define the starting permutation
    ind_start = 1 + (c - 1) * core_perm_n;
    
    % Define the last permutation
    ind_end = ind_start + core_perm_n - 1;
    
    % Generate a permutation indicator that signals each core which
    % permutations it needs to process
    core_perm_cell{c} =  ind_start:ind_end;
    
end
% Core Loop

% Create a matrix containing the maximum cluster t sums for each
% permutation
ctsmax_perm = zeros(2, core_perm_n);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Prepare a function for the permutation procedure
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Function that allows the parallel caluclation of permutations
function [ctsmax_perm] = chanfun(c, core_perm_cell, data_perm, shiftmat, ...
    vpn, chans, times, freqs, nghbrhd, conn, dropchan, ctsmax_perm)
    
    % Permutation Loop
    for p = 1:length(core_perm_cell{c})
    

        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
        %% Permute averages
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

    
        % Duplicate the original data matrix
        power_shift = data_perm;
        
        % Extract the permutation indicator
        ind_perm = core_perm_cell{c}(p);
                
        % Swap averages around
        power_shift(shiftmat(:,ind_perm),:,:) = ...
            power_shift(shiftmat(:,ind_perm),[2,1],:);
            
            
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
        %% Perform a baseline normalization
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
        
        
        % Reshape data
        power_shift = reshape(power_shift, [length(vpn), 2, ...
            length(chans), length(times), length(freqs)]);
            
        % Calculate the baseline
        bl = repmat(mean(power_shift,4), [1,1,1,length(times),1]);

        % Baseline the data
        power_bl = 100 .* ((power_shift - bl) ./ bl);
            
        % Calculate the condition difference
        power_diff = squeeze(power_bl(:,1,:,:,:) - power_bl(:,2,:,:,:));
        
        % Reshape data
        power_diff = reshape(power_diff, [length(vpn), ...
            length(chans) * length(times) * length(freqs)]);
        
        
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
        %% Calculate t values
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
        
        
        % Calculate t values against zero
        power_tval = ttest0(power_diff);
        
        % Reshape data
        power_tval = reshape(power_tval, ...
            [length(chans), length(times), length(freqs)]);
                
        
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
        %% Organize channel data topographically
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %        
        
        
        % Topographical x loop
        for x = 1:size(nghbrhd,1)

            % Topographical y loop
            for y = 1:size(nghbrhd,2)

                % Check whether the topographical spot is a channel
                if nghbrhd(x,y) > 0

                    % Copy data at the right topographical spot
                    topo_tval(x,y,:,:) = ...
                        squeeze(power_tval(nghbrhd(x,y),:,:));

                end
                % Check whether the topographical spot is a channel

            end
            % Topographical Y Loop

        end
        % Topographical X Loop
        
        
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
        %% Peform the clustering procedure
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
        
        
        % Determine significant condition differences
        topo_sig = topo_tval;
        topo_sig(topo_tval > 2.1314) = 1;
        topo_sig(topo_tval < -2.1314) = -1;


        % Sign Loop
        for s = 1:2
        
            % Extract sign values of the same direction
            if s == 1
                topo_sign = topo_sig == -1;
            else
                topo_sign = topo_sig == 1;
            end
            % Extract sign values of the same direction

            % Cluster neighbors
            clustered_mat = bwlabeln(topo_sign, conn);

            % Correct for channels that occured double
            for i = 1:size(dropchan,2)
                clustered_mat(dropchan(1,i),dropchan(2,i),:,:) = 0;
            end
            
            % Predefine a vector containg all summed up cluster t
            % statistics for this permutation
            cts_perm = zeros(1,max(clustered_mat(:)));
            
            % Empty Statement
            if isempty(cts_perm) == 0
            
                % Cluster Loop
                for l = 1:max(clustered_mat(:))
                    
                    % Create a cluster indicator
                    ind = clustered_mat == l;
                    
                    % Calculate the cluster t sums
                    cts_perm(1,l) = sum(topo_tval(ind));
                    
                end
                % Cluster Loop
                
                % Extract the largest or smallest cluster t sum
                if s == 1
                    ctsmax_perm(s,p) = min(cts_perm);
                else
                    ctsmax_perm(s,p) = max(cts_perm);
                end
                % Extract the largest or smallest cluster t sum
                
            else
            
                % If there are no significant values the cluster t sum
                % is set to 0
                ctsmax_perm(s,p) = 0;
                
            end
            % Empty Statement
            
        end
        % Sign Loop
                
        % Display progress
        if c == 1
            disp(['Permutation ', num2str(p), '/', ...
                num2str(size(ctsmax_perm,2)), ...
                ' has been calculated.']);
        end
        % Display progress        
        
    end
    % Permutation Loop
    
endfunction
% Function that allows the parallel caluclation of thresholding parameters
% by distributing the channel data over cores

% Define a channel variable whose expressions are distributed over the
% cores
parachanfun = @(c)chanfun(c, core_perm_cell, data_perm, shiftmat, ...
    vpn, chans, times, freqs, nghbrhd, conn, dropchan, ctsmax_perm);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Perform the thresholding
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

    
% Start the permutations
[ctsmax_cell] = pararrayfun(cores_n, parachanfun, cores_ind, ...
    'UniformOutput', false);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Transform the cell into a matrix structure
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Cell Loop
for i = 1:size(ctsmax_cell,2)

    if i == 1
        ctsmax_mat = ctsmax_cell{1,i};
    else
        ctsmax_mat = [ctsmax_mat, ctsmax_cell{1,i}];
    end
    
end
% Cell Loop


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Save the cluster t sums
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Save the thresolding data
save('-mat7-binary', [wd, 'results/', studlab, '/tfa/power/', ...
    'global/maris_og/spec_bl/power_pm_ctsmax.mat'], ...
    'ctsmax_mat');


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Prepare the clustering procedure for the true data
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Calculate the condition difference
power_diff = squeeze(power_tm(:,1,:,:,:) - power_tm(:,2,:,:,:));

% Reshape data
power_diff = reshape(power_diff, [length(vpn), ...
    length(chans) * length(times) * length(freqs)]);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Calculate t values
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Calculate t values against zero
power_tval = ttest0(power_diff);

% Reshape data
power_tval = reshape(power_tval, ...
    [length(chans), length(times), length(freqs)]);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Organize channel data topographically
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %        


% Topographical x loop
for x = 1:size(nghbrhd,1)

    % Topographical y loop
    for y = 1:size(nghbrhd,2)

        % Check whether the topographical spot is a channel
        if nghbrhd(x,y) > 0

            % Copy data at the right topographical spot
            topo_tval(x,y,:,:) = ...
                squeeze(power_tval(nghbrhd(x,y),:,:));

        end
        % Check whether the topographical spot is a channel

    end
    % Topographical Y Loop

end
% Topographical X Loop


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Peform the clustering procedure
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Determine significant condition differences
topo_sig = topo_tval;
topo_sig(topo_tval > 2.1314) = 1;
topo_sig(topo_tval < -2.1314) = -1;


% Sign Loop
for s = 1:2

    % Extract sign values of the same direction
    if s == 1
        topo_sign = topo_sig == -1;
    else
        topo_sign = topo_sig == 1;
    end
    % Extract sign values of the same direction

    % Cluster neighbors
    clustered_mat = bwlabeln(topo_sign, conn);

    % Correct for channels that occured double
    for i = 1:size(dropchan,2)
        clustered_mat(dropchan(1,i),dropchan(2,i),:,:) = 0;
    end

    % Predefine a matrix in which cluster indicators will be replaced
    % by their cluster t sum (cts)
    if s == 1
        cts_topo = zeros(size(clustered_mat));
    end

    % Cluster Loop
    for l = 1:max(clustered_mat(:))

        % Create a cluster indicator
        ind = clustered_mat == l;

        % Calculate the cluster t sums
        cts_topo(ind) = sum(topo_tval(ind));

    end
    % Cluster Loop

end
% Sign Loop


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
cts_mat = zeros(length(chans), length(times), length(freqs));

% Topographical x loop
for x = 1:size(nghbrhd_red,1)

    % Topographical y loop
    for y = 1:size(nghbrhd_red,2)
    
        % Only copy data if the spatial position corresponds to a channel
        if nghbrhd_red(x,y) > 0
        
            % Copy data from a topographical to a normal channel
            % structure
            cts_mat(nghbrhd_red(x,y),:,:) = squeeze(cts_topo(x,y,:,:));
            
        end
        % Only copy data if the spatial position corresponds to a channel
    
    end
    % Topographical y loop
    
end
% Topographical x loop


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% Save the cluster t sums
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Save the thresolding data
save('-mat7-binary', [wd, 'results/', studlab, '/tfa/power/', ...
    'global/maris_og/spec_bl/power_tm_cts.mat'], ...
    'cts_mat', 'chans', 'times', 'freqs');
    
    
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%% End: Conduct inference statistics according to Maris & Oostenveld
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

