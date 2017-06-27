function [] = runERFPlots(subject)

addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
ft_defaults
megpath = '/mnt/homes/home024/chernandez/meg_data/surprise/preprocessed/Data/ERF/';  % path of preprocessed MEG data

cfg = data.cfg;
cfg.vartrllength = 1;
avg =ft_timelockanalysis(cfg, data);
cfgtopo.xlim = [0.3 0.5];
cfg.colorbar = 'yes';
ft_topoplotER(cfg, avg);




%%%%%%%%%%%%%%%%%%%%%%%%%
load('CTF275_helmet.mat');
lay.outline = lay.outline([1 3:end]); % remove outer bound

cfgtopo                     = [];
cfgtopo.marker              = 'off';
cfgtopo.layout              = lay;
cfgtopo.comment             = 'no';
cfgtopo.highlight           = 'on';
cfgtopo.highlightsymbol     = '.';
cfgtopo.highlightsize       = 3;
cfgtopo.highlightchannel    = data.label;%chans(c).names;
%cfgtopo.colormap            = cmap;
%cfgtopo.zlim                = [0 6e-14];%conditions(n).crange; % use dB instead of %
cfgtopo.xlim = [0.3 0.5]; %TODO calculate the interval -100 ms before gocue
cfgtopo.renderer            = 'painters';

%topodata.dimord     = 'chan_freq_time';
topodata.dimord     = 'subj_chan__time';
cfgtopo.style           = 'straight_imsat';
ft_topoplotER(cfgtopo, topodata);



end