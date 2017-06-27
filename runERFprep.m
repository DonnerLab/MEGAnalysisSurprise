function data = runERFprep(subject,session,recording)
  % Low-pass filter ERF data (40Hz), downsample (100Hz), apply planar
  % gradient transformation, and save

% ==================================================================
% SPECIFY PATHS AND GET SUBJECT-SPECIFIC FILES
% ==================================================================
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
ft_defaults

megpath = '/mnt/homes/home024/chernandez/meg_data/surprise/preprocessed/Data/';  % path of preprocessed MEG data
behavpath = ['/mnt/homes/home024/pmurphy/Surprise_accumulation/Data/',subject,filesep];  % path of behavioural data
savepath = '/mnt/homes/home024/chernandez/meg_data/surprise/preprocessed/Data/ERF/';


% addpath D:\Experiments\Surprise_accumulation\Analysis\Behaviour
% addpath D:\Experiments\Surprise_accumulation\Analysis\Gen_fun
% addpath 'C:\Program Files\MATLAB\fieldtrip-20160221'  % tell Matlab where FieldTrip is
% ft_defaults
% 
% megpath = 'D:\Experiments\Surprise_accumulation\Analysis\MEG\Data\';
% behavpath = ['D:\Experiments\Surprise_accumulation\Data\',subject,filesep];
% savepath = 'D:\Experiments\Surprise_accumulation\Analysis\MEG\Data\ERF\';

subjfile = dir([megpath,filesep,subject,'-',session,'*',recording,'.mat']);  % pull current meg filename

% ==================================================================
% LOAD CURRENT MEG/BEHAVIOURAL DATA & DERIVE COMPUTATIONAL VARIABLES
% ==================================================================
load([megpath,subjfile.name])  % load meg data

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
% LOW-PASS FILTER ERFs
% ==================================================================
% low-pass filter settings take to match de Lange et al. (2010), J Neurosci
cfg             = [];
cfg.lpfilter    = 'yes'; % request low-pass filter
cfg.lpfreq      = 40; % low-pass cutoff in hz
cfg.lpfiltord   = 6;  % filter order
cfg.lpfilttype  = 'but';  % filter type
cfg.lpfiltdir   = 'twopass';
data            = ft_preprocessing(cfg, data);

% ==================================================================
% DOWNSAMPLE
% ==================================================================
fsnew = 100;  % desired sampling rate

% adjust all sample references in event matrix to fit new sampling rate
data.trialinfo(:,[3:4 7:10]) = round(data.trialinfo(:,[3:4 7:10])./data.fsample.*fsnew);

% resample
cfg             = [];
cfg.resamplefs  = fsnew;
cfg.detrend     = 'no';
data = ft_resampledata(cfg, data);

% ==================================================================
% PLANAR GRADIENT TRANSFORMATION
% ==================================================================
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

% compute planar gradient magnitude over both directions
data = ft_combineplanar([], data);

% ==================================================================
% SAVE
% ==================================================================
data.mdlvars.LLR = LLR_full;
data.mdlvars.LPR = LPR_full;
data.mdlvars.surprise = surprise_full;
data.mdlvars.deltaL = deltaL_full;
data.mdlvars.choices = choices_full;

save([savepath,subject,'-',session,'-',recording,'_ERF.mat'],'data')
    
end