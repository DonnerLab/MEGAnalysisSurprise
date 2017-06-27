function [] = runERFanalysis(subject)
  % Applies two baselining operations (pre-trial and peri-sample) to MEG
  % ERF data and regresses these data (channel*time) onto model-based
  % variables of interest

basewin1 = [-0.1 0];  % baseline relative to pre-mask onset (s)
basewin2 = [0 0.05];  % baseline relative to each sample onset (s)

% ==================================================================
% SPECIFY PATHS AND GET SUBJECT-SPECIFIC FILES
% ==================================================================
% addpath /mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
ft_defaults
megpath = '/mnt/homes/home024/chernandez/meg_data/surprise/preprocessed/Data/ERF/';  % path of preprocessed MEG data

% addpath 'C:\Program Files\MATLAB\fieldtrip-20160221'  % tell Matlab where FieldTrip is
% ft_defaults
% megpath = 'D:\Experiments\Surprise_accumulation\Analysis\MEG\Data\ERF\';

subjfiles = dir([megpath,'/',subject,'*.mat']);  % pull all meg filenames for this subject

smp_data1=[]; smp_data2=[]; go_data1=[]; go_data2=[];
LLR_full=[]; LPR_full=[]; surprise_full=[]; deltaL_full=[]; choices_full=[]; sess_full=[];
for f = 1:length(subjfiles)

    load([megpath,subjfiles(f).name])  % load meg data
    sess = str2double(subjfiles(f).name(5));
    
    % ==================================================================
    % APPLY BASELINE I (PRE-MASK ONSET)
    % ==================================================================
    fprintf('Applying baseline I...\n')
    for t = length(data.trial)
        stbase = mean(data.trial{t}(:,data.time{t}>=basewin1(1) & data.time{t}<=basewin1(2)),2);
        data.trial{t} = data.trial{t}-repmat(stbase,1,size(data.trial{t},2));
    end
    
    % ==================================================================
    % APPLY BASELINE II (PER-SAMPLE ONSET) & CONCATENATE SEGMENTS
    % ==================================================================
    fprintf('Applying baseline II and concatenating sample-wise data segments...\n')
    onsets = 0.4:0.4:0.4*13;  % vector of all sample onset times relative to pre-mask
    if f==1, slen = length(find(data.time{1}>=onsets(1) & data.time{1}<onsets(2))); end  % getting fixed sample length so there's no small timing conflicts across trials
    for t = 1:length(data.trial)
        tc = size(smp_data1,1)+1;
        for s = 1:length(onsets)-1
            stsmp = find(data.time{t}>=onsets(s),1,'first');
            smp_data1(tc,:,:,s) = data.trial{t}(:,stsmp:stsmp+slen-1); % store sample-segmented data from baseline regime I
            
            stbase = mean(data.trial{t}(:,data.time{t}>=onsets(s)+basewin2(1) & data.time{t}<=onsets(s)+basewin2(2)),2);
            smp_data2(tc,:,:,s) = ...   % store sample-segmented data from baseline regime II
                data.trial{t}(:,stsmp:stsmp+slen-1)-repmat(stbase,1,length(stsmp:(stsmp+slen-1)));  % applying sample-wise baseline
        end
        go_data1(end+1,:,:) = data.trial{t}(:,(end-data.fsample*0.6):end);  % store pre-go cue data from baseline regime I
        
        stbase = mean(data.trial{t}(:,data.time{t}>=onsets(end)+basewin2(1) & data.time{t}<=onsets(end)+basewin2(2)),2);
        go_data2(end+1,:,:) = data.trial{t}(:,(end-data.fsample*0.6):end)-repmat(stbase,1,length(data.time{t}((end-data.fsample*0.6):end)));  % store pre-go cue data from baseline regime II
    end
    
    % ==================================================================
    % CONCATENATE MODEL-BASED VARIABLES
    % ==================================================================
    LLR_full = [LLR_full; data.mdlvars.LLR];
    LPR_full = [LPR_full; data.mdlvars.LPR];
    surprise_full = [surprise_full; data.mdlvars.surprise];
    deltaL_full = [deltaL_full; data.mdlvars.deltaL];
    choices_full = [choices_full; data.mdlvars.choices];
    sess_full = [sess_full; ones(length(data.mdlvars.choices),1).*sess];
    
    % ==================================================================
    % STORE FT STRUCTURES FOR LATER PLOTTING
    % ==================================================================
    grad = data.grad;
    cfg = data.cfg;
end

% ==================================================================
% RUN SINGLE-TRIAL REGRESSIONS
% ==================================================================
% Baseline approach 1
Tprior1=[]; Tposterior1=[]; Tchoice1=[];
fprintf('Running regressions of ERF onto prior belief (trial-wise baseline)...\n')
for s = 1:size(smp_data1,4)-1  % looping through samples
    for c = 1:size(smp_data1,2)  % looping through channels
        for t = 1:size(smp_data1,3)  % looping through time-points
            % Regressing ERF data (ERF:s+1) onto prior belief (L:s)
            m = regstats(smp_data1(:,c,t,s+1),LPR_full(:,s),'linear',{'tstat'});
            Tprior1(c,s,t) = m.tstat.t(2);
        end
    end
end
fprintf('Running regressions of ERF onto posterior belief (trial-wise baseline)...\n')
for s = 1:size(smp_data1,4)  % looping through samples
    for c = 1:size(smp_data1,2)  % looping through channels
        for t = 1:size(smp_data1,3)  % looping through time-points
            % Regressing ERF data (ERF:s) onto posterior belief (L:s)
            m = regstats(smp_data1(:,c,t,s),LPR_full(:,s),'linear',{'tstat'});
            Tposterior1(c,s,t) = m.tstat.t(2);
        end
    end
end
fprintf('Running regressions of ERF onto choice (trial-wise baseline)...\n')
for c = 1:size(go_data1,2)  % looping through channels
    for t = 1:size(go_data1,3)  % looping through time-points
        % Regressing ERF data onto choice
        m = regstats(go_data1(:,c,t),choices_full,'linear',{'tstat'});
        Tchoice1(c,t) = m.tstat.t(2);
    end
end

% Baseline approach 2
Tprior2=[]; Tposterior2=[]; Tsurprise2=[]; TdeltaL2=[]; Tchoice2=[];
fprintf('Running regressions of ERF onto prior belief, surprise and belief change (sample-wise baseline)...\n')
for s = 1:size(smp_data2,4)-1  % looping through samples
    for c = 1:size(smp_data2,2)  % looping through channels
        for t = 1:size(smp_data2,3)  % looping through time-points
            m = regstats(smp_data2(:,c,t,s+1),LPR_full(:,s),'linear',{'tstat'});
            Tprior2(c,s,t) = m.tstat.t(2);
            
            m = regstats(smp_data2(:,c,t,s+1),surprise_full(:,s+1),'linear',{'tstat'});
            Tsurprise2(c,s,t) = m.tstat.t(2);
            
            m = regstats(smp_data2(:,c,t,s+1),abs(deltaL_full(:,s+1)),'linear',{'tstat'});
            TdeltaL2(c,s,t) = m.tstat.t(2);
        end
    end
end
fprintf('Running regressions of ERF onto posterior belief (sample-wise baseline)...\n')
for s = 1:size(smp_data2,4)  % looping through samples
    for c = 1:size(smp_data2,2)  % looping through channels
        for t = 1:size(smp_data2,3)  % looping through time-points
            m = regstats(smp_data2(:,c,t,s),LPR_full(:,s),'linear',{'tstat'});
            Tposterior2(c,s,t) = m.tstat.t(2);
        end
    end
end
fprintf('Running regressions of ERF onto choice (sample-wise baseline)...\n')
for c = 1:size(go_data2,2)  % looping through channels
    for t = 1:size(go_data2,3)  % looping through time-points
        % Regressing ERF data onto choice
        m = regstats(go_data2(:,c,t),choices_full,'linear',{'tstat'});
        Tchoice2(c,t) = m.tstat.t(2);
    end
end

% ==================================================================
% SAVE RESULTS
% ==================================================================
save([megpath,subject,'_regression_output.mat'],'Tprior1','Tposterior1','Tchoice1','Tprior2','Tsurprise2','TdeltaL2','Tposterior2','Tchoice2','grad','cfg')





