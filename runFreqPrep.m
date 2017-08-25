function [] = runFreqPrep(subject,session,recording)
  % Apply planar gradient transformation and calculate TF representations of
  % MEG power for each of two gradiometers per sensor/trial, using a single
  % Hanning taper for freq range 3-35Hz (window length: 400 ms, step size:
  % 50 ms, freq resolution: 2.5 Hz, bin size: 1 Hz). TF representations are calculated both aligned to trial onset and  After TFR calculation,
  % power estimates from the two planar gradiometers for each sensor are
  % combined by taking their sum.

% ==================================================================
% SPECIFY PATHS AND GET SUBJECT-SPECIFIC FILES
% ==================================================================
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
ft_defaults

megpath = '/mnt/homes/home024/chernandez/meg_data/surprise/preprocessed/Data/';  % path of preprocessed MEG data
behavpath = ['/mnt/homes/home024/pmurphy/Surprise_accumulation/Data/',subject,filesep];  % path of behavioural data
savepath = '/mnt/homes/home024/chernandez/meg_data/surprise/preprocessed/Data/TF/';

subjfile = dir([megpath,filesep,subject,'-',session,'*',recording,'.mat']);  % pull current meg filename

% ==================================================================
% LOAD CURRENT MEG/BEHAVIOURAL DATA & DERIVE COMPUTATIONAL VARIABLES
% ==================================================================
fprintf('\nLoading meg file: %s...\n',subjfile.name)
load([megpath,subjfile.name])  % load meg data

fprintf('\nLoading and processing behavioural data...\n')
stimIn_full=[]; LLR_full=[]; LPR_full=[]; surprise_full=[]; deltaL_full=[]; choices_full=[];
for b = unique(data.trialinfo(:,1))'  % looping through each block within this meg dataset
    load([behavpath,'S',session,filesep,'Behaviour',filesep,subject,'_',session,'_',num2str(b),'.mat'])
    load([behavpath,'S',session,filesep,'Sample_seqs',filesep,subject,'_',session,'_',num2str(b),'.mat'])
    
    Behav = Behav(unique(data.trialinfo(data.trialinfo(:,1)==b,2)),:);  % dumping any trials not contained in meg dataset
    stimIn = stimIn(unique(data.trialinfo(data.trialinfo(:,1)==b,2)),:);
    
    % Converting sample and choice values to appropriate signs for choice regressions
    stimIn = round(stimIn.*-1);
    choices = Behav(:,2)-1;
    
    % Convert stimulus values to LLRs & calculate sample-wise surprise
    LLRin = log(normpdf(stimIn,gen.mu(2)*-1,gen.sigma(2))./normpdf(stimIn,gen.mu(1)*-1,gen.sigma(1)));
    LPR=[]; surprise=[]; deltaL=[];
    for t = 1:length(choices)
        [LPR(t,1:size(LLRin,2)),surprise(t,1:size(LLRin,2))] = accGlaze(LLRin(t,:),gen.H,0,'DY');
        [~,deltaL(t,1:size(LLRin,2))] = accGlaze(LLRin(t,:),gen.H,0,'absL');
    end
    
    % Collating useable single trials
    stimIn_full = [stimIn_full; stimIn];        % sample sequences 
    LLR_full = [LLR_full; LLRin];               % sample evidence strength
    LPR_full = [LPR_full; LPR];                 % evolving belief
    surprise_full = [surprise_full; surprise];  % surprise
    deltaL_full = [deltaL_full; deltaL];        % change in belief
    choices_full = [choices_full; choices];     % trial-by-trial choices
end

% ==================================================================
% PULL ONLY FULL-LENGTH TRIALS - consider skipping this to maximize trial counts
% ==================================================================
fprintf('Keeping only full-length trials...\n')
assert(length(choices_full)==length(data.trial),'ERROR: Trial counts in MEG/behaviour are unequal')

ts=[];  % useable trials based on behavioural data
nsamps=[];
for t = 1:length(choices_full)
    nsamps(t,1) = length(find(~isnan(stimIn_full(t,:))));
    if sum(isnan(stimIn_full(t,:)))==0 && choices_full(t)<2, ts(end+1) = t; end
end
assert(isempty(find((data.trialinfo(:,end)-nsamps)~=0, 1)),'ERROR: Mismatch in MEG/behaviour number of samples per trial')

LLR_full = LLR_full(ts,:);
LPR_full = LPR_full(ts,:);
surprise_full = surprise_full(ts,:);
deltaL_full = deltaL_full(ts,:);
choices_full = choices_full(ts,:);

cfg             = [];
cfg.trials      = ts;
data = ft_selectdata(cfg, data);

% ==================================================================
% PLANAR GRADIENT TRANSFORMATION
% ==================================================================
fprintf('\nRunning planar gradient transformation...\n')

% define neighbours based on CTF template
cfg                 = [];
cfg.method          = 'template';
cfg.layout          = 'CTF275';
neighbours          = ft_prepare_neighbours(cfg);

% compute planar gradiometers for MEG sensors
cfg                 = [];
cfg.feedback        = 'no';
cfg.method          = 'template';
cfg.planarmethod    = 'sincos';
cfg.channel         = 'MEG';
cfg.neighbours      = neighbours;
data                = ft_megplanar(cfg, data);

% ==================================================================
% TIME-FREQUENCY DECOMPOSITION
% ==================================================================
fprintf('\nRunning time-frequency decomposition...\n')

cfg                 = [];
cfg.output          = 'pow';
cfg.channel         = 'MEG';
cfg.method          = 'mtmconvol';   % specifying multi-taper method
cfg.taper           = 'hanning';     % with Hanning taper
cfg.keeptrials      = 'yes';
cfg.keeptapers      = 'no';
cfg.precision       = 'single'; % saves disk space
%cfg.feedback        = 'none'; % improve readability of logfiles

% make nice timebins at each 50 ms, will include 0 point of pre-mask onset; last time point will correspond to go cue time for trial with shortest dot-to-go interval
mintime = data.time{1}(1);
minpos = find(data.trialinfo(:,3)==min(data.trialinfo(:,3))); maxtime = data.time{minpos(1)}(end);
toi = floor(mintime) : 0.05 : ceil(maxtime);
toi(toi < mintime) = []; toi(toi > maxtime) = [];

cfg.toi             = toi; % all time within each locking
%cfg.pad             = 4; % pad to a fixed number of seconds before TFR

cfg.foi             = 3:1:35;   % frequencies of interest
cfg.t_ftimwin       = ones(1, length(cfg.foi)) .* 0.4;
% cfg.t_ftimwin       = 3*(1./cfg.foi);  % adaptive time window of 3 cycles per estimated frequency
cfg.tapsmofrq       = [];

freq                 = ft_freqanalysis(cfg, data);
% assert(isequal(freq.freq, cfg.foi), '! spectral estimation returned different foi, double-check settings !');

% ==================================================================
% COMBINE PLANAR GRADIOMETERS
% ==================================================================
freq = ft_combineplanar([], freq);

% ==================================================================
% SAVE
% ==================================================================
freq.mdlvars.LLR = LLR_full;
freq.mdlvars.LPR = LPR_full;
freq.mdlvars.surprise = surprise_full;
freq.mdlvars.deltaL = deltaL_full;
freq.mdlvars.choices = choices_full;

save([savepath,subject,'-',session,'-',recording,'_TF.mat'],'freq')
    
end