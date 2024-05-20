% D409V/WT medial septum MEA analysis 
% Author: Millie Sander
% Part 1: data processing 

% filter to get only high quality cells: 
p_bursts_bl_hq = p_bursts_bl(ind);
p_bursts_ss_hq = p_bursts_ss(ind);

% p(burst) 
function [p_bursts] = isi_outlier_loopBG(units_sel, phase_starts, phase_ends)
%Summary of this function goes here
%   Calculate the waveform symmetry index for all units' average waveform
%DEPENDENCIES:
% wv_symmetry.mlx (function written by JB for undergrad2020
%%
%Preallocate
p_bursts = NaN(length(units_sel), length(phase_ends));
%Calculate burst indices 
for i = 1:length(units_sel)
    for ii = 1:length(phase_ends)
    p_bursts(i, ii) = isi_outlier(diff(units_sel(i).tstamp(units_sel(i).tstamp < phase_ends(ii) & units_sel(i).tstamp > phase_starts(ii))));
    end
end
end

[p_bursts_bl] = isi_outlier_loopBG(units_sel, info.timecourse(1).timeWin(1, 13)-180, info.timecourse(1).timeWin(1, 13))

[p_bursts_ss] = isi_outlier_loopBG(units_sel, info.timecourse(1).timeWin(2, 13)-180, info.timecourse(1).timeWin(2, 13))

% raster plot 
function [h1] = rasterplot_one(filename, phase_end, phasewin, ch)
% This function will make raster plots for all units in this experiment. 
% INPUTS:
%       filename = '210517Aa.mat'
%       phasewin = time (seconds) frame you want to plot the rasterplots
%       n        = neuron number (specify which neuron you want to plot)
% OUTPUTS:
%       h1       = handle for figure
%
%%
clc
%close all
load(filename);
%%
phase_start = phase_end - phasewin;
%%
%Set parameters for subsequent rasterplots
%Define the segment of time of the baseline you want to look at.
final_ts = cell(1,1);
nspikes = nan(1,1);
% Create a new timestamp (t) vector variable with elements that falls
% between baseline_start and baseline_end.
final_ts{1,1} = units_sel(ch).tstamp(units_sel(ch).tstamp < phase_end & units_sel(ch).tstamp > phase_start);
nspikes = numel(final_ts);
%
%%
%Raster plot
%Each line (subplot) in the figure is a putative neuron. The vertical lines represent spikes. 
scrsz = get(0,'ScreenSize');
h1 = figure('Position',scrsz)
for i=1:length(final_ts)
    rasterSpikes = final_ts{i};
    plot(length(final_ts))
    set(gca,'xlim',[phase_start phase_end])
    for ii = 1:numel(rasterSpikes)
        line([rasterSpikes(ii) rasterSpikes(ii)] ,[0 10] , 'Color','k');
    end
    %set(gca,'ytick',[])
    %set(gca,'yticklabel',[])
    % or
    axis off 
end
xlabel('Time (s)','fontsize',20, 'visible','on');
ylabel(strcat('Neuron', " " ,num2str(ch)),'fontsize',20, 'rotation',90,'horizontalalignment','left', 'verticalalignment','middle', 'visible','on');
end

% for entire 210722b_MS experiment, neuron 1 
[h1_one_peak] = rasterplot_one('210722b_MS', 5445, 5445, 1)
% for 3 mins of baseline, then 10 seconds around peak, and 3 mins of ss
[h1_one_peak] = rasterplot_one('210722b_MS', 1800, 180, 1)
[h1_one_peak] = rasterplot_one('210722b_MS', pooled_sel.histHz.t(timepointnew(1))-5, 10, 1)
[h1_one_peak] = rasterplot_one('210722b_MS', 3600, 180, 1)

% isi 
function [isi, isi_mean] = isi_MS(units_sel, phase_starts, phase_ends)
%Summary of this function goes here
%   Calculate the waveform symmetry index for all units' average waveform
%DEPENDENCIES:
% wv_symmetry.mlx (function written by JB for undergrad2020
%%
%Preallocate
isi = cell(length(units_sel), 1);
%Calculate burst indices 
for i = 1:length(units_sel)
    for ii = 1:length(phase_ends)
    isi{i, 1} = diff(units_sel(i).tstamp(units_sel(i).tstamp < phase_ends(ii) & units_sel(i).tstamp > phase_starts(ii)));
    end
    isi_mean(i) = mean(isi{i, 1});
end

end

% timeWin(condition, recording)
% so (1, 2) is baseline cut off for second experiment, 
% and (2, 5) is drug cut off for 5th recording 

[isi_bl, isi_mean_bl] = isi_MS(new_units_sel, info.timecourse(1).timeWin(1, 13)-180, info.timecourse(1).timeWin(1, 13))
[isi_ss, isi_mean_ss] = isi_MS(new_units_sel, info.timecourse(1).timeWin(2, 13)-180, info.timecourse(1).timeWin(2, 13))

% [isi_test, isi_mean_test] = isi_MS(units_sel, phasestarts, phaseends);

[firing_rates_sel, timepoint] = mea_peak_firing_rate_ind(pooled_sel,S.drugs,firing_rates_sel);

timepointnew= timepoint(ind);
for n = 1: size(new_units_sel, 2)
phasestarts (n) = pooled_sel.histHz.t(timepointnew(n))-15+info.timecourse(1).timeWin (1,2);
phaseends(n) = pooled_sel.histHz.t(timepointnew(n))+15+info.timecourse(1).timeWin (1,2);
[isi_peak, isi_mean_peak] = isi_MS(new_units_sel, phasestarts(n), phaseends(n));
end



 
