
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Plot ERP results for the original Maris & Oostenveld approach
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %


% Clear the workspace
clear all;
close all hidden;


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Load and define variables

% Laden der stabilen Variablen
wd = 'D:\Experimental\TempInt_EEG\';

% Berücksichtigung weiterer Funktionen
addpath([wd, 'functions\altmany-export_fig-9d97e2c']);

% Open eeglab
eeglab;

% Define the participant variable
vpn = [2,3,6,7,12,14,15,16,19,22,26,27,28,30,33];

% Load the cluster results of the true data
clust_tm = getfield(load([wd, '\results\experimental\erp\global\', ...
    'maris_og\pos_tm_cts.mat'], 'cts_mat'), 'cts_mat');

% Load the cluster results of the bootstrapped data
clust_pm = getfield(load([wd, '\results\experimental\erp\global\', ...
    'maris_og\pos_pm_ctsmax.mat'], 'ctsmax_mat'), 'ctsmax_mat');

% Participant Loop
for v = 1:length(vpn)
    
    % Load raw erp data
    erp_tm(v,:,:,:) = mean(getfield(load([wd, '\results\experimental\', ...
        'erp\raw\v', num2str(vpn(v)), '_epoched_bp.1-50.mat'], ...
        'erp_ep'), 'erp_ep'),4);
    
end
% Participant Loop

% Load channel information
chans = getfield(load([wd, '\results\experimental\', ...
        'erp\raw\v', num2str(vpn(1)), '_epoched_bp.1-50.mat'], ...
        'chans'), 'chans');

% Load time information
times = getfield(load([wd, '\results\experimental\', ...
        'erp\raw\v', num2str(vpn(1)), '_epoched_bp.1-50.mat'], ...
        'times'), 'times');

% Load topographical information
chantopo_num = getfield(load([wd, '\results\', ...
    'clustering_fundamentals.mat'], 'chantopo_num'), 'chantopo_num');

% Load channel information
chanlocs = getfield(load([wd, '\results\', ...
    'chanlocs.mat'], 'chanloc'), 'chanloc');


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Determine significant clusters

% Sort the bootstrapped values
clust_pms = sort(abs(clust_pm(:)));

% Define the critical value
critval = clust_pms(round(.95 * length(clust_pms)));

% Set insignificant clusters to zero
clust_tm(clust_tm > (-1 * critval) & clust_tm < critval) = 0;


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Plot the topography of each cluster

% Calculate the average participant difference
erp_avg = squeeze(mean(erp_tm(:,1,:,:) - erp_tm(:,2,:,:),1));

% Define cluster parameters
ct = {[200,300],[300,400]};

% Define cluster labels
clustlab = {'t200-300', 't300-400'};

% Cluster Loop
for c = 1:length(ct)
    
    % Define cluster indicator variables
    times_ind = times >= ct{c}(1) & times <= ct{c}(2);
    
    % Create a plotting matrix
    erp_plot = mean(erp_avg(:,times_ind),2);
    
    % Reduce time
    times_indred = times_ind(times >= 0);
    
    % Determine significant channels
    sigchansind = ...
        find(sum(clust_tm(:,times_indred) < 0,2) > 0);
    
    % Open figure
    figure('unit','centimeters','position',[2,2,13,13]);
    
    % Plot the topography
    topoplot(erp_plot, chanlocs, 'conv', 'on', ...
        'emarker2', {sigchansind, '.', 'w', 25, 1});

    % Change colormap
    colormap(jet);

    % Change layout of the colorbar
    cb = colorbar('vert', 'FontName', 'Arial', ...
        'Fontsize', 12, 'FontWeight', 'normal', ...
        'location', 'southoutside');
    set(cb, 'YDir' ,'normal');

    % Scale the colorbar
    caxis([-1.25,1.25]);
    set(cb, 'ylim', [-1.25,1.25], 'ytick', -1:.5:1, 'yticklabel', {}, ...
        'FontName', 'Arial', ...
        'FontSize', 12, 'FontWeight', 'normal');

    % Rescale image
    set(gca, 'unit', 'centimeters', 'position', [.5,.5,12,12]);

    % Rescale colorbar
    set(cb, 'visible', 'off');
    
    % Save figure
    export_fig([wd, 'images\erp\experimental\maris_og\global\', ...
        'manuscript\topoplot_', clustlab{c}], ...
        '-bmp', '-nocrop', '-r500');
    
    % Close figure
    close;
    
end
% Cluster Loop


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Plot the grandaverage erp curves


%% Prepare data to calculate the confidence intervals

% Restructure matrix
data_ConVp = permute(erp_tm, [2,1,3,4]);

% Average data over conditions but not participants
data_Vp = repmat(mean(data_ConVp,1),size(data_ConVp,1),1,1,1);

% Additionally average over participants
data = repmat(mean(data_Vp,2),1,size(data_Vp,2),1,1);

% Calculate data free from between-subject variance (Cousineau, 2005)
data_wsv = data_ConVp - data_Vp + data;


%% Calculate within-subject confidence intervals

% Calculate the participant average
mean_wsv = squeeze(mean(data_wsv,2));

% Calculate the participant standard deviation
std_wsv = squeeze(std(data_wsv, 0, 2));

% Calculate the number of participants
n_wsv = size(data_wsv,2);

% Determine the critical t value for calculating the cofidence intervals
tCI = tinv(.975, size(data_wsv,2)-1);

% Calculate the upper and lower confidence interval boundaries
cis_wsv = mean_wsv - tCI .* ((std_wsv) ./ sqrt(n_wsv));
cil_wsv = mean_wsv + tCI .* ((std_wsv) ./ sqrt(n_wsv));


%% Average significant channel

% Define the time window of interest
reltime = times >= 200 & times <= 400;

% Reduce time
reltime_red = reltime(times >= 0);

% Detect significant channels in the time window of interest
sigchan_ind = ...
    squeeze(sum(clust_tm(:,reltime_red) < 0, 2)) > 51;

% Only use the average significiant channel for plotting the later ERP
erp_mean = squeeze(mean(mean_wsv(:,sigchan_ind,:),2));
erp_cis = squeeze(mean(cis_wsv(:,sigchan_ind,:),2));
erp_cil = squeeze(mean(cil_wsv(:,sigchan_ind,:),2));


%% Plot the grand average ERP condition courses
    
% Open figure
figure('Units', 'centimeters', 'position', [2,2,18,11], 'color', 'w');

% Define colors (1 = simultaneity (red), 2 = separation (blue))
colormap = {[1,0,0], [0,1,0]};
linestyle = {'-','-.'};
colormap_fill = {[1,.5,.5], [.5,1,.5]};

% Condition loop
for b = 1:size(erp_mean,1)

    % Extract data to plot
    dplotMean = squeeze(erp_mean(b,:));
    dplotCILow = squeeze(erp_cis(b,:));
    dplotCIHigh = squeeze(erp_cil(b,:));

    % Plot data
    plot(times, dplotMean, 'color', colormap{b}, ...
        'LineStyle', linestyle{b}, 'LineWidth', 2);
    hold on;
    f = fill([times, fliplr(times)], ...
        [dplotCILow, fliplr(dplotCIHigh)], colormap_fill{b});
    set(f, 'facealpha', .4, 'edgecolor', 'none', ...
        'HandleVisibility', 'off');

end
% Condition loop

% Adapt axes layout
set(gca, 'ylim', [-3,6.5], 'ytick', -2:2:6, 'yticklabel', {}, ...
    'FontName', 'Times New Roman', 'Fontsize', 10, 'FontWeight', 'normal');
set(gca, 'xlim', [-200,900], 'xtick', -100:100:800, 'xticklabel', {}, ...
    'FontName', 'Times New Roman', 'Fontsize', 10, 'FontWeight', 'normal');

% Add line with respect to the target presentation and 0 on the y-axis
line([0,0],[-3,6.5], 'color', [0,0,0], 'LineStyle', '-', ...
    'HandleVisibility', 'off');
line([min(times),max(times)],[0,0], 'color', [0,0,0], ...
    'LineStyle', '-', 'HandleVisibility', 'off');

% Highlight significant time interval
f = fill([200, 450, 450, 200], [-3, -3, 6.5, 6.5], [.75,.75,.75]);
set(f, 'facealpha', .25, 'edgecolor', 'none', 'HandleVisibility', 'off');

% Rescale image
set(gca, 'unit', 'centimeters', 'position', [.5,.5,17,10]);

% Save the figure as a .bmp-file
export_fig([wd, 'images\erp\experimental\maris_og\global\', ...
    'manuscript\grandaverage_sigchans'], ...
    '-bmp', '-nocrop', '-r500');

% Close figure
close;


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Plot ERP results for the original Maris & Oostenveld approach
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

