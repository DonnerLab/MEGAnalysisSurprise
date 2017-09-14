function [] = runFreqAnalysis_LIonly(subject)
  % Applies two baselining operations (pre-trial and peri-sample) to MEG
  % ERF data and regresses these data (channel*time) onto model-based
  % variables of interest

basewin = [-0.3 -0.2];  % baseline window relative to pre-mask onset (s)
smpwin = [0 0.8];  % window for sample-wise analyses
trlwin = [-0.3 5.8];  % window for full-trial-wise analyses
respwin = [-0.2 0];  % window relative to freq data end for response prep analysis

Lchans = {'MLC21','MLC13','MLC22','MLC41','MLF63','MLC14','MLC23','MLC31','MLF54','MLF64','MLC15','MLC24','MLF45','MLF55','MLF65','MLC16','MLF46','MLF56','MLF66'}; % left motor/premotor channels
Rchans = {'MRC21','MRC13','MRC22','MRC41','MRF63','MRC14','MRC23','MRC31','MRF54','MRF64','MRC15','MRC24','MRF45','MRF55','MRF65','MRC16','MRF46','MRF56','MRF66'}; % right motor/premotor channels

% ==================================================================
% SPECIFY PATHS AND GET SUBJECT-SPECIFIC FILES
% ==================================================================
addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
ft_defaults
megpath = '/mnt/homes/home024/chernandez/meg_data/surprise/preprocessed/Data/TF/';  % path of preprocessed MEG data

subjfiles = dir([megpath,subject,'-*.mat']);  % pull all meg filenames for this subject

smp_LI=[]; zsmp_LI=[]; trl_LI=[]; go_data=[];
LLR_full=[]; LPR_full=[]; surprise_full=[]; deltaL_full=[]; choices_full=[]; sess_full=[]; samp_full=[]; priors_full=[];

% ==================================================================
% TRIAL-WISE ANALYSIS
% ==================================================================
fprintf('Beginning trial-wise analysis...\n')
for f = 1:length(subjfiles)

    load([megpath,subjfiles(f).name])  % load meg data
    sess = str2double(subjfiles(f).name(5));
    
    % ==================================================================
    % STORE FT STRUCTURES FOR LATER PLOTTING
    % ==================================================================
    if f==1,
        grad = freq.grad;
        cfg = freq.cfg;
        
        Lchans = find(ismember(grad.label,Lchans));  % converting channel strings into channel indices
        Rchans = find(ismember(grad.label,Rchans));  % converting channel strings into channel indices
    end
    
    % ==================================================================
    % APPLY BASELINES
    % ==================================================================
    fprintf('Applying baseline for %s...\n',subjfiles(f).name)
    
    freq.time = round(freq.time,2);  % rounding time vector to nearest 2nd decimal - otherwise slight inaccuracies can lead to bad timing later
    %freq.time = round_ndec(freq.time,2); % ch: since round(x,n) is not supported for version 2013a
    stbase = nanmean(10.*log10(freq.powspctrm(:, :, :, freq.time>=basewin(1) & freq.time<=basewin(2))), 4);  % extract dB-transformed pre-onset baselines per trial, channel & freq
    for t = 1:size(stbase,1)   % looping through trials
        for c = 1:size(stbase,2)  % looping through channels
            freq.powspctrm(t,c,:,:) = 10.*log10(squeeze(freq.powspctrm(t,c,:,:)))-repmat(squeeze(stbase(t,c,:)),1,length(freq.time));  % applying baselines to dB-transformed data
        end
    end
    
    % ==================================================================
    % PULL DESIRED SEGMENTS OF DATA & CALCULATE LATERALIZATION INDEX
    % ==================================================================
    fprintf('Concatenating trial- & sample-wise data segments...\n')
    if f==1,
        trltimes = freq.time(freq.time>=trlwin(1) & freq.time<=trlwin(2));
    end
    for t = 1:size(freq.powspctrm,1)
        trl_LI(end+1,:,:) = squeeze(nanmean(freq.powspctrm(t,Rchans,:,freq.time>=trlwin(1) & freq.time<=trlwin(2)))) - ...
            squeeze(nanmean(freq.powspctrm(t,Lchans,:,freq.time>=trlwin(1) & freq.time<=trlwin(2))));  % store full-trial LI data - trials*freqs*times
        go_data(end+1,:,:) = nanmean(freq.powspctrm(t,:,:,freq.time>=(max(freq.time)+respwin(1)) & freq.time<=(max(freq.time)+respwin(2))),4);  % store time-averaged post-sequence, pre-go-cue data
    end
    
    % ==================================================================
    % CONCATENATE MODEL-BASED VARIABLES
    % ==================================================================
    LLR_full = [LLR_full; freq.mdlvars.LLR];
    LPR_full = [LPR_full; freq.mdlvars.LPR];
    surprise_full = [surprise_full; freq.mdlvars.surprise];
    deltaL_full = [deltaL_full; freq.mdlvars.deltaL];
    choices_full = [choices_full; freq.mdlvars.choices];
    sess_full = [sess_full; ones(length(freq.mdlvars.choices),1).*sess];
    samps=1:1:12;
    for i=1:length(freq.mdlvars.choices)
        for s=1:length(samps)
            samp_full = [samp_full;s];
            if s~=1; priors_full = [priors_full;s];end
        end
    end
end

% ==================================================================
% MAKE CATEGORICAL SESSION REGRESSORS
% ==================================================================
sessions = unique(sess_full);
sess_r = zeros(length(sess_full),length(sessions)-1);
for s = 1:length(sessions)-1
    sess_r(sess_full==sessions(s+1),s) = ones(length(find(sess_full==sessions(s+1))),1);
end
% categorical variable for stacked samples analysis
samps = unique(samp_full);
samp_r = zeros(length(samp_full),length(samps)-1);
for s = 1:length(samps)-1
    samp_r(samp_full==samps(s+1),s) = ones(length(find(samp_full==samps(s+1))),1);
end
% categorical variable for stacked samples analysis (prior and surprise analysis)
samps = unique(priors_full);
priors_r = zeros(length(priors_full),length(samps)-1);
for s = 1:length(samps)-1
    priors_r(priors_full==samps(s+1),s) = ones(length(find(priors_full==samps(s+1))),1);
end
size(priors_r)
% ==================================================================
% COMPUTE TRIAL-AVERAGED RESPONSES
% ==================================================================
resp_avg = squeeze(mean(go_data(choices_full==1,:,:),1))-squeeze(mean(go_data(choices_full==0,:,:),1));  % chan*freq matrix of average right minus left responses in post-sequence, pre-go period

% ==================================================================
% RUN SINGLE-TRIAL REGRESSIONS
% ==================================================================
% Full-trial regressions, lateralization index
TposteriorF_LI=[]; TchoiceF_LI=[];
fprintf('Running regressions of full-trial lateralization index (time*freq) onto final choice and final posterior belief...\n')
for f = 1:size(trl_LI,2)  % looping through freqs
    for t = 1:size(trl_LI,3)  % looping through time-points
        m = regstats(trl_LI(:,f,t),[choices_full sess_r],'linear',{'tstat'});  % final choice
        TchoiceF_LI(f,t) = m.tstat.t(2);
        m = regstats(trl_LI(:,f,t),[LPR_full(:,end) sess_r],'linear',{'tstat'});  % posterior belief
        TposteriorF_LI(f,t) = m.tstat.t(2);
    end
end

% ==================================================================
% PULL FINAL-SAMPLE SINGLE-TRIAL LATERALIZATION, Z-SCORED WITHIN SESSION
% ==================================================================
timewin = [4.7 4.9];
finalLI = nan(size(trl_LI,1),size(trl_LI,2));
for f = 1:size(trl_LI,2);
    for s = sessions';
        finalLI(sess_full==s,f) = squeeze(zscore(mean(trl_LI(sess_full==s,f,trltimes>=timewin(1) & trltimes<=timewin(2)),3)));
    end
end
final_prior = LPR_full(:,end-1);

% ==================================================================
% SAVE RESULTS AND CLEAN UP
% ==================================================================
freqs = freq.freq;
save([megpath,subject,'_trialwise_output_LIonly.mat'],'trltimes','freqs','grad','cfg','resp_avg','TchoiceF_LI','TposteriorF_LI','finalLI','final_prior')

clear trl_data go_data trl_LI resp_avg TchoiceF_LI TposteriorF_LI finalLI final_prior


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ==================================================================
% SAMPLE-WISE ANALYSIS
% ==================================================================
fprintf('Beginning sample-wise analysis...\n')
for f = 1:length(subjfiles)

    load([megpath,subjfiles(f).name])  % load meg data
    
    % ==================================================================
    % APPLY BASELINES
    % ==================================================================
    fprintf('Applying baseline for %s...\n',subjfiles(f).name)
    
    freq.time = round(freq.time,2);  % rounding time vector to nearest 2nd decimal - otherwise slight inaccuracies can lead to bad timing later
    
    stbase = nanmean(10.*log10(freq.powspctrm(:, :, :, freq.time>=basewin(1) & freq.time<=basewin(2))), 4);  % extract dB-transformed pre-onset baselines per trial, channel & freq
    for t = 1:size(stbase,1)   % looping through trials
        for c = 1:size(stbase,2)  % looping through channels
            freq.powspctrm(t,c,:,:) = 10.*log10(squeeze(freq.powspctrm(t,c,:,:)))-repmat(squeeze(stbase(t,c,:)),1,length(freq.time));  % applying baselines to dB-transformed data
        end
    end
    
    % ==================================================================
    % PULL DESIRED SEGMENTS OF DATA
    % ==================================================================
    fprintf('Concatenating trial- & sample-wise data segments...\n')
    onsets = 0.4:0.4:0.4*12;  % vector of all sample onset times relative to pre-mask
    if f==1,
        smptimes = freq.time(freq.time>=smpwin(1) & freq.time<=smpwin(2));  % getting vector of sample times relative to dot onset
    end
    for t = 1:size(freq.powspctrm,1)
        tc = size(smp_LI,1)+1;
        for s = 1:length(onsets)
            stsmp = find(freq.time>=(onsets(s)+smpwin(1)),1,'first');
            zerosmp = find(freq.time<onsets(s),1,'last');
            smp_LI(tc,:,:,s) = squeeze(nanmean(freq.powspctrm(t,Rchans,:,stsmp:stsmp+length(smptimes)-1),2)) - ... % store sample-segmented data (trial*freq*time*sample)
                squeeze(nanmean(freq.powspctrm(t,Lchans,:,stsmp:stsmp+length(smptimes)-1),2)); % right minus left so right choices will have positive values (same as LLR signing)
            zsmp_LI(tc,:,s) = squeeze(nanmean(freq.powspctrm(t,Rchans,:,zerosmp),2)) - ... % store power from time point preceding each sample onset (trial*freq*sample)
                squeeze(nanmean(freq.powspctrm(t,Lchans,:,zerosmp),2)); % right minus left so right choices will have positive values (same as LLR signing)
            
        end
    end
end
clear freq

% ==================================================================
% RUN SINGLE-TRIAL REGRESSIONS
% ==================================================================
% Sample-wise regressions, lateralization index
TpriorS_LI=[]; TpriorWllrS_LI=[]; TpriorllrS_LI=[]; TllrXsurpriseS_LI=[];
fprintf('Running regressions of sample-wise lateralization index (time*freq) onto signed prior and LLR*surprise...\n')
for s = 1:size(smp_LI,4)-1  % looping through samples
    for f = 1:size(smp_LI,2)  % looping through freqs
        for t = 1:size(smp_LI,3)  % looping through time-points
            try
                m = regstats(smp_LI(:,f,t,s+1),[LPR_full(:,s) sess_r],'linear',{'tstat'});  % prior
                TpriorS_LI(f,t,s) = m.tstat.t(2);
                m = regstats(smp_LI(:,f,t,s+1),[LPR_full(:,s) nanzscore(LLR_full(:,s+1)) nanzscore(LLR_full(:,s+1)).*nanzscore(surprise_full(:,s+1)) sess_r],'linear',{'tstat'});  % LLR and LLR*surprise
                TpriorWllrS_LI(f,t,s) = m.tstat.t(2);
                TpriorllrS_LI(f,t,s) = m.tstat.t(3);
                TllrXsurpriseS_LI(f,t,s) = m.tstat.t(4);
            catch ME
                TpriorS_LI(f,t,s) = nan;
                TpriorWllrS_LI(f,t,s) = nan;
                TpriorllrS_LI(f,t,s) = nan;
                TllrXsurpriseS_LI(f,t,s) = nan;
            end
        end
    end
end
TposteriorS_LI=[]; TllrS_LI=[];
fprintf('Running regressions of sample-wise lateralization index (time*freq) onto posterior and LLR...\n')
for s = 1:size(smp_LI,4)  % looping through samples
    for f = 1:size(smp_LI,2)  % looping through freqs
        for t = 1:size(smp_LI,3)  % looping through time-points
            try
                m = regstats(smp_LI(:,f,t,s),[LPR_full(:,s) sess_r],'linear',{'tstat'});  % posterior
                TposteriorS_LI(f,t,s) = m.tstat.t(2);
                m = regstats(smp_LI(:,f,t,s),[LLR_full(:,s) sess_r],'linear',{'tstat'});  % LLR
                TllrS_LI(f,t,s) = m.tstat.t(2);
            catch ME
                TposteriorS_LI(f,t,s) = nan;
                TllrS_LI(f,t,s) = nan;
            end
        end
    end
end

TprecWllrS_LI=[]; TprecllrS_LI=[]; TprecllrXsurpriseS_LI=[];
fprintf('Running regressions of sample-wise lateralization index (time*freq) onto precedent timestamp power and LLR*surprise...\n')
for s = 1:size(smp_LI,4)-1  % looping through samples
    for f = 1:size(smp_LI,2)  % looping through freqs
        for t = 1:size(smp_LI,3)  % looping through time-points
            try
                m = regstats(smp_LI(:,f,t,s+1),[zsmp_LI(:,f,s+1) nanzscore(LLR_full(:,s+1)) nanzscore(LLR_full(:,s+1)).*nanzscore(surprise_full(:,s+1)) sess_r],'linear',{'tstat'});  % LLR and LLR*surprise
                TprecWllrS_LI(f,t,s) = m.tstat.t(2);
                TprecllrS_LI(f,t,s) = m.tstat.t(3);
                TprecllrXsurpriseS_LI(f,t,s) = m.tstat.t(4);
            catch ME
                ME.message
                TprecWllrS_LI(f,t,s) = nan;
                TprecllrS_LI(f,t,s) = nan;
                TprecllrXsurpriseS_LI(f,t,s) = nan;
            end
        end
    end
end
% ==================================================================
% SAVE RESULTS
% ==================================================================
save([megpath,subject,'_samplewise_output_LIonly.mat'],'smptimes','freqs','grad','cfg',...
    'TpriorS_LI','TpriorWllrS_LI','TpriorllrS_LI','TllrXsurpriseS_LI','TposteriorS_LI','TllrS_LI',...
    'TprecWllrS_LI','TprecllrS_LI','TprecllrXsurpriseS_LI')


% ==================================================================
% RUN SAMPLE-WISE REGRESSIONS, 'STACKED' APPROACH
% ==================================================================
% Sample-wise regressions, lateralization index
LPR_fsmp = []; LPR_fprior = []; LLR_fsmp = []; smp_LI_fsmp = []; zsmp_LI_fsmp = []; surprise_fsmp = [];

% stack samples for prior belief and LLR*surprise analysis
for s = 1:size(smp_LI,4)
    if s==1; LPR_fprior = [LPR_fprior; LPR_full(:,s)]; end
    if s~=1
        smp_LI_fsmp = [smp_LI_fsmp; smp_LI(:,:,:,s)]; 
        zsmp_LI_fsmp = [zsmp_LI_fsmp; zsmp_LI(:,:,s)]; 
        LLR_fsmp = [LLR_fsmp; LLR_full(:,s)];
        surprise_fsmp = [surprise_fsmp; surprise_full(:,s)];    
    end
    if s~=12 && s~=1; LPR_fprior = [LPR_fprior; LPR_full(:,s)]; end
end

% adjust categorical variables for prior and surprise sample analysis
fsp_sess_r = [];
for i=1:length(sess_r)
    for s=1:12
        if s~=1; fsp_sess_r=[fsp_sess_r;sess_r(i,:)]; end
    end
end


TpriorS_LI=[]; TpriorWllrS_LI=[]; TpriorllrS_LI=[]; TllrXsurpriseS_LI=[];
fprintf('Running regressions of sample-wise lateralization index (time*freq) onto signed prior and LLR*surprise, stacked approach...\n')
for f = 1:size(smp_LI_fsmp,2)  % looping through freqs
    for t = 1:size(smp_LI_fsmp,3)  % looping through time-points
        try
            m = regstats(smp_LI_fsmp(:,f,t),[LPR_fprior fsp_sess_r priors_r],'linear',{'tstat'});  % prior
            TpriorS_LI(f,t) = m.tstat.t(2);
            m = regstats(smp_LI_fsmp(:,f,t),[LPR_fprior nanzscore(LLR_fsmp) nanzscore(LLR_fsmp).*nanzscore(surprise_fsmp) fsp_sess_r priors_r],'linear',{'tstat'});  % LLR and LLR*surprise
            TpriorWllrS_LI(f,t) = m.tstat.t(2);
            TpriorllrS_LI(f,t) = m.tstat.t(3);
            TllrXsurpriseS_LI(f,t) = m.tstat.t(4);
        catch ME
            TpriorS_LI(f,t) = nan;
            TpriorWllrS_LI(f,t) = nan;
            TpriorllrS_LI(f,t) = nan;
            TllrXsurpriseS_LI(f,t) = nan;
        end
    end
end

% Sample-wise regressions, lateralization index, using the power of the
% precedent sample, stacked approach
TprecWllrS_LI=[]; TprecllrS_LI=[]; TprecllrXsurpriseS_LI=[];
fprintf('Running regressions of sample-wise lateralization index (time*freq) onto precedent timestamp power and LLR*surprise, stacked approach...\n')
for f = 1:size(smp_LI,2)  % looping through freqs
    for t = 1:size(smp_LI,3)  % looping through time-points
        try
            m = regstats(smp_LI_fsmp(:,f,t),[zsmp_LI_fsmp(:,f) nanzscore(LLR_fsmp) nanzscore(LLR_fsmp).*nanzscore(surprise_fsmp) fsp_sess_r priors_r],'linear',{'tstat'});  % LLR and LLR*surprise
            TprecWllrS_LI(f,t) = m.tstat.t(2);
            TprecllrS_LI(f,t) = m.tstat.t(3);
            TprecllrXsurpriseS_LI(f,t) = m.tstat.t(4);
        catch ME
            ME.message
            TprecWllrS_LI(f,t) = nan;
            TprecllrS_LI(f,t) = nan;
            TprecllrXsurpriseS_LI(f,t) = nan;
        end
    end
end

% stack sample variables for posterior belief
smp_LI_fsmp=[]; LPR_fsmp=[]; LLR_fsmp=[];
for s = 1:size(smp_LI,4)
    smp_LI_fsmp = [smp_LI_fsmp; smp_LI(:,:,:,s)];
    LPR_fsmp = [LPR_fsmp; LPR_full(:,s)];
    LLR_fsmp = [LLR_fsmp; LLR_full(:,s)];
end
% adjust categorical variable session
fs_sess_r = [];
for i=1:length(sess_r)
    for j=1:12
        fs_sess_r=[fs_sess_r;sess_r(i,:)];
    end
end

TposteriorS_LI=[]; TllrS_LI=[];
fprintf('Running regressions of sample-wise lateralization index (time*freq) onto posterior and LLR, stacked approach...\n')
for f = 1:size(smp_LI_fsmp,2)  % looping through freqs
    for t = 1:size(smp_LI_fsmp,3)  % looping through time-points
        try
            m = regstats(smp_LI_fsmp(:,f,t),[LPR_fsmp fs_sess_r samp_r],'linear',{'tstat'});  % posterior
            TposteriorS_LI(f,t) = m.tstat.t(2);
            m = regstats(smp_LI_fsmp(:,f,t),[LLR_fsmp fs_sess_r samp_r],'linear',{'tstat'});  % LLR
            TllrS_LI(f,t) = m.tstat.t(2);
        catch ME
            TposteriorS_LI(f,t,s) = nan;
            TllrS_LI(f,t) = nan;
        end
    end
end

save([megpath,subject,'_samplewise_stacksam_LIonly.mat'],'smptimes','freqs','grad','cfg',...
    'TpriorS_LI','TpriorWllrS_LI','TpriorllrS_LI','TllrXsurpriseS_LI','TposteriorS_LI','TllrS_LI',...
    'TprecWllrS_LI', 'TprecllrS_LI', 'TprecllrXsurpriseS_LI')
end




