%{

March 15th, 2024

Original script my by Ilhan Bok + TDT.
Modified by Ido Haber for Psychoplastogens project @ Hai's Lab UW Madison


Running this script will:

1. Create a .csv with processed data for analysis
2. Figures
    - Histogram per dataset
    - Raster per dataset
    - Histogram per pair (baseline - intervention)
    - Histogram averaged across groups.
    - Summary figures for freq. & active channels. 

%}

% 1. Housekeeping
close all; clear; clc;

% 2. Set initial parameters
SDKPATH = '/Volumes/ahai/TDT/TDTMatlabSDK';

%ketanserin path:
%BASEPATH = '/Volumes/ahai/Shared/Psychoplastogens Project/MEA data/00_GoodDatasetsForAnalysis/Ketanserin_11012023';

%DOI path:
BASEPATH = '/Volumes/ahai/Shared/Psychoplastogens Project/MEA data/00_GoodDatasetsForAnalysis/NewBatch_09140223';

desktopPath = '/Users/idohaber/Desktop/Final_DOI'; % Set this to your output directory
addpath(genpath(SDKPATH));

% Parameters for data processing
STORE = 'Wav1';
CHANNELS = 1:64; % channels we want to process for raster
BP_FILTER = [300 2500]; % Band-pass filter parameters
f0 = 60; % Notch filter center frequency
fs = 48000; % Sampling frequency
fn = fs/2;
freqRatio = f0/fn;
[bw, a] = iirnotch(freqRatio, freqRatio/35);

% Ketanserin sets:
% datasets = {...
%     'IdoControl-231101-144046_#1', ...
%     'IdoKetanserin-231101-133725_#1', ...
%     'IdoControl-231101-145330_#2', ...
%     'IdoKetanserin-231101-134946_#2', ...
%     'IdoControl-231101-150512_#3', ...
%     'IdoKetanserin-231101-140237_#3'};

% DOI sets:
datasets = {...
    'IdoControl-230914-130200_#1', ...
    'IdoDOI-230914-142502_#1', ...
    'IdoControl-230914-131548_#2', ...
    'IdoDOI-230914-143740_#2', ...
    'IdoControl-230914-132855_#3', ...
    'IdoDOI-230914-144945_#3', ...
    'IdoControl-230914-154022_#4', ...
    'IdoDOI-230914-161838_#4', ...
    'IdoControl-230914-155318_#5', ...
    'IdoDOI-230914-163140_#5', ...
    'IdoControl-230914-160601_#6', ...
    'IdoDOI-230914-164512_#6'}; 


% Initialize metrics storage
activeChannelsPerDataset = zeros(length(datasets), 1);
totalSpikesPerDataset = zeros(length(datasets), 1);
firingFrequencyPerDataset = zeros(length(datasets), 1);
spikeCountsPerChannelPerDataset = zeros(length(CHANNELS), length(datasets));

% Define a table to store dataset parameters
paramsTable = table('Size', [length(datasets) 5], ...
    'VariableTypes', {'string', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'DatasetName', 'FiringFrequency Spike/Sec', 'TotalSpikes', '# of ActiveChannels', 'RecordLengthSecs'});

% Loop through each dataset
for i = 1:length(datasets)
    ENDPATH = datasets{i};
    BLOCKPATH = fullfile(BASEPATH, ENDPATH);

    datasetSpikes = 0;
    datasetActiveChannels = 0;
    recordLengthSecs = 0;
    spikeCountsPerChannel = zeros(length(CHANNELS), 1);

    rasterData = cell(length(CHANNELS), 1);
    for ch = CHANNELS
        data = TDTbin2mat(BLOCKPATH, 'STORE', STORE, 'CHANNEL', ch);
        dataFiltered = TDTdigitalfilter(data, STORE, BP_FILTER, 'ORDER', 5);
        dataFiltered.streams.(STORE).data = filtfilt(bw, a, dataFiltered.streams.(STORE).data);
        
        spikes = TDTthresh(dataFiltered, STORE, 'MODE', 'auto', 'POLARITY', -1, 'STD', 6.5, 'TAU', 5);
        rasterData{ch} = spikes.snips.Snip.ts;

        if ch == CHANNELS(1)
            fs = data.streams.(STORE).fs;
            nSamples = length(data.streams.(STORE).data);
            recordLengthSecs = nSamples / fs;
        end

        if ~isempty(rasterData{ch})
            spikeCount = length(rasterData{ch});
            datasetSpikes = datasetSpikes + spikeCount;
            spikeCountsPerChannel(ch) = spikeCount;
            datasetActiveChannels = datasetActiveChannels + 1;
        end
    end

    % Store metrics
    activeChannelsPerDataset(i) = datasetActiveChannels;
    totalSpikesPerDataset(i) = datasetSpikes;
    firingFrequencyPerDataset(i) = datasetSpikes / recordLengthSecs;
    spikeCountsPerChannelPerDataset(:, i) = spikeCountsPerChannel;

    % Store metrics in the table
    paramsTable.DatasetName(i) = ENDPATH;
    paramsTable.FiringFrequency(i) = firingFrequencyPerDataset(i);
    paramsTable.TotalSpikes(i) = totalSpikesPerDataset(i);
    paramsTable.ActiveChannels(i) = activeChannelsPerDataset(i);
    paramsTable.RecordLengthSecs(i) = recordLengthSecs;

    % Plotting Raster Plot
    figure; hold on;
    for ch = 1:length(CHANNELS)
        spikeTimes = rasterData{ch};
        for spike = 1:length(spikeTimes)
            plot([spikeTimes(spike), spikeTimes(spike)], [ch-0.4, ch+0.4], 'k');
        end
    end
    title(sprintf('%s\nActive Channels = %d, Total Spikes = %d, Firing Frequency = %.2f spikes/s', ENDPATH, datasetActiveChannels, datasetSpikes, firingFrequencyPerDataset(i)));
    xlabel('Time (s)'); ylabel('Channel'); axis tight;
    saveas(gcf, fullfile(desktopPath, sprintf('raster_%s.png', ENDPATH)), 'png');
    close;

    % Spike Count per Channel Plot
    figure;
    bar(1:length(CHANNELS), spikeCountsPerChannel, 'FaceColor', 'b');
    title(sprintf('Spike Count per Channel for %s', ENDPATH));
    xlabel('Channel Number'); ylabel('Spike Count');
    grid on; axis tight;
    saveas(gcf, fullfile(desktopPath, sprintf('spike_count_per_channel_%s.png', ENDPATH)), 'png');
    close;
end

% After the loop, save the table to a .csv file
writetable(paramsTable, fullfile(desktopPath, 'dataset_parameters.csv'));


%% Compute differences and ratios for data storage and Visualization 


% Initialize arrays to store differences and ratios for each pair
numPairs = length(datasets) / 2;
spikesPerChannelDiff = zeros(length(CHANNELS), numPairs);
firingFreqRatio = zeros(numPairs, 1);
activeChannelsDiff = zeros(numPairs, 1);

% Compute differences and ratios for each pair
for i = 1:numPairs
    % Indices for odd (1, 3, 5, ...) and even (2, 4, 6, ...) datasets
    oddIdx = 2*i - 1;
    evenIdx = 2*i;
    
    % Compute differences in spikes per channel
    spikesPerChannelDiff(:, i) = spikeCountsPerChannelPerDataset(:, oddIdx) - spikeCountsPerChannelPerDataset(:, evenIdx);

    % Compute ratios in firing frequency
    firingFreqRatio(i) = firingFrequencyPerDataset(oddIdx) / firingFrequencyPerDataset(evenIdx);

    % Compute differences in active channels
    activeChannelsDiff(i) = activeChannelsPerDataset(oddIdx) - activeChannelsPerDataset(evenIdx);
end

% Compute the average spike count difference
averageSpikeCountPerChannelOdd = mean(spikeCountsPerChannelPerDataset(:, 1:2:end), 2);
averageSpikeCountPerChannelEven = mean(spikeCountsPerChannelPerDataset(:, 2:2:end), 2);
averageSpikeCountDiff = averageSpikeCountPerChannelOdd - averageSpikeCountPerChannelEven;


writetable(paramsTable, fullfile(desktopPath, 'dataset_parameters.csv'));



%% Plotting:

% Plotting Figure 1 with color-coded bars
figure;
subplot(3,1,1);
hold on;
for i = 1:length(CHANNELS)
    if averageSpikeCountDiff(i) > 0
        bar(i, averageSpikeCountDiff(i), 'FaceColor', "#A2142F");
    else
        bar(i, averageSpikeCountDiff(i), 'FaceColor', "#0072BD");
    end
end
hold off;
title('Average Δ in Spike Counts per Channel (Baseline - DOI + Ket)');
xlabel('Channel Number'); ylabel('Average Δ in Spike Count');

% Rotate x-axis labels and adjust font size
xticks(1:length(averageSpikeCountDiff)); % Set x-axis ticks to show integers
xtickangle(90); % Rotate x-axis labels to 90 degrees
set(gca, 'FontSize', 8); % Set font size to 8 (or any other suitable size)


% b. Ratio in firing frequency
subplot(3, 1, 2);
hold on;
for i = 1:length(firingFreqRatio)
    if firingFreqRatio(i) > 1
        bar(i, firingFreqRatio(i), 'FaceColor', "#A2142F");
    else
        bar(i, firingFreqRatio(i), 'FaceColor', "#0072BD");
    end
end
hold off;
title('Ratio in Firing Frequency: Baseline / DOI + Ket');
xlabel('Dataset Pair'); ylabel('Firing Freq. Ratio');
set(gca, 'XTick', 1:numPairs, 'XTickLabel', {'1/2', '3/4', '5/6'});

% c. Difference in active channels
subplot(3, 1, 3);
hold on;
for i = 1:length(activeChannelsDiff)
    if activeChannelsDiff(i) > 0
        bar(i, activeChannelsDiff(i), 'FaceColor', "#A2142F");
    else
        bar(i, activeChannelsDiff(i), 'FaceColor', "#0072BD");
    end
end
hold off;
title('Δ in Active Channels: Baseline - DOI + Ket');
xlabel('Dataset Pair'); ylabel('Δ in Active Channels');
set(gca, 'XTick', 1:numPairs, 'XTickLabel', {'1-2', '3-4', '5-6'});

saveas(gcf, fullfile(desktopPath, 'figure1_comparisons.png'), 'png');

%% Figures for Differences between baseline -> intervention

% Generate figures for the difference in spike count per channel for each pair with color coding
for pairIndex = 1:numPairs
    figure;
    hold on; % Allow multiple bars to be plotted on the same figure
    
    % Iterate over each channel to plot and color each bar based on its value
    for channelIndex = 1:length(CHANNELS)
        barColor = '#0072BD'; % Default to blue for negative differences
        if spikesPerChannelDiff(channelIndex, pairIndex) > 0
            barColor = '#A2142F'; % Change to red for positive differences
        end
        bar(channelIndex, spikesPerChannelDiff(channelIndex, pairIndex), 'FaceColor', barColor);
    end
    
    hold off;
    title(sprintf('Difference in Spike Counts per Channel (Pair %d: %s - %s)', ...
        pairIndex, datasets{2*pairIndex-1}, datasets{2*pairIndex}));
    xlabel('Channel Number');
    ylabel('Difference in Spike Count');
    
    % Adjustments for readability
    xticks(1:length(CHANNELS)); % Show all channel numbers
    xtickangle(90); % Rotate labels for better readability
    set(gca, 'FontSize', 8); % Adjust font size if necessary
    
    % Save each figure with a unique name based on the pair index
    saveas(gcf, fullfile(desktopPath, sprintf('spike_count_diff_pair_%d_color_coded.png', pairIndex)), 'png');
end
