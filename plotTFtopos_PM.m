clear, close all

addpath /mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221 
addpath /mnt/homes/home024/chernandez/matlab/cbrewer
addpath /mnt/homes/home024/chernandez/matlab/MEGAnalysisSurprise
% addpath D:\Experiments\Surprise_accumulation\Analysis\Gen_fun
ft_defaults
loadpath = '/mnt/homes/home024/chernandez/meg_data/surprise/preprocessed/Data/TF/';  % path of preprocessed MEG data

% load some stuff for plotting
load('CTF275_helmet.mat');  % load layout
lay.outline = lay.outline([1 3:end]); % remove outer bound

[T,cmap] = evalc(['cbrewer(''div'', ''RdBu'', 64)']);  % create colourmap (and suppressing warning messages using evalc) - Caro, you need to download the cbrewer set of functions for this to work...
cmap = flipud(cmap); % red = increase                  % ... get those here: https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab

%subjects
allsubj = {'DHB','EXF','TFD'};
allsubj = {'DCB' 'DHB' 'ECB' 'EMB' 'EXF' 'EXG' 'GSB' 'HBC' 'JTB' 'KSV' 'NIF' 'OMF' 'PDP' 'QNV' 'TFD' 'TNB' 'TSJ' };

% frequency bands of interest
alphabnd = [8 14];
betabnd = [15 28];
allbnd = [8 28];

% create cfg structure for topo plotting
cfgtopo                     = [];
cfgtopo.marker              = 'off';
cfgtopo.layout              = lay;
cfgtopo.comment             = 'no';
%cfgtopo.highlight           = 'on';
%cfgtopo.highlightsymbol     = '.';
%cfgtopo.highlightsize       = 3;
%cfgtopo.highlightchannel    = chans(c).names;
cfgtopo.colormap            = cmap;
cfgtopo.renderer            = 'painters';
cfgtopo.style               = 'straight_imsat';

topodata.dimord     = 'chan_time';
topodata.time       = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- PLOT CHANNEL-AVERAGED TF RESPONSE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Load and process data
% trialTF=[];
% for s = 1:length(allsubj)
%     load([loadpath,allsubj{s},'_trialwise_output.mat'])
%     trialTF(:,:,s) = mean(trl_avg(:,freqs>=5,:),1);
% end
% trialTF(:,:,end+1) = mean(trialTF,3);
% 
% % Plot
% onsets = 0.4:0.4:(0.4*12);
% figure, 
% for s = 1:length(allsubj)+1
%     subplot(length(allsubj)+1,1,s), hold on
%     imagesc(trltimes,freqs(freqs>=5),trialTF(:,:,s),[-max(max(abs(trialTF(:,:,s)))) max(max(abs(trialTF(:,:,s))))]),
%     line([0 0],[min(freqs(freqs>=5)) max(freqs)],'LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
%     line([onsets; onsets],[repmat(min(freqs(freqs>=5)),1,length(onsets)); repmat(max(freqs(freqs>=5)),1,length(onsets))],'LineStyle','--','Color',[0.3 0.3 0.3])
%     set(gca,'YDir','normal'), xlim([min(trltimes) max(trltimes)]); ylim([min(freqs(freqs>=5)) max(freqs)])
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- 1)  PLOT LATERALIZATION TOPOGRAPHIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and process data
alphatopo=[]; betatopo=[]; alltopo=[];
for s = 1:length(allsubj)
    load([loadpath,allsubj{s},'_trialwise_output_LIonly.mat'])
    alphatopo(:,s) = mean(resp_avg(:,freqs>=alphabnd(1) & freqs<=alphabnd(2)),2);
    betatopo(:,s) = mean(resp_avg(:,freqs>=betabnd(1) & freqs<=betabnd(2)),2);
    alltopo(:,s) = mean(resp_avg(:,freqs>=allbnd(1) & freqs<=allbnd(2)),2);
end
alphatopo(:,end+1) = mean(alphatopo,2);
betatopo(:,end+1) = mean(betatopo,2);
alltopo(:,end+1) = mean(alltopo,2);

% Processing topo stuff
freqs = round(freqs);

grad.chanpos = grad.chanpos(1:274,:);
grad.chanori = grad.chanori(1:274,:);
grad.chantype = grad.chantype(1:274);
grad.chanunit = grad.chanunit(1:274);
grad.label = grad.label(1:274);

topodata.label      = grad.label;
topodata.grad       = grad;

% Plot
spos = [1 4 7 10];
snfig = [5 9 13 17];
%spos = 1:3:51;
sp=0;
figure, 
for s = 1:length(allsubj)+1
    if ~isempty(find(snfig==s))
        figure, sp=1;
    else
        sp=sp+1;
    end
    %subplot(length(allsubj)+1,3,spos(sp)), hold on
    subplot(4,3,spos(sp)), hold on
    topodata.avg = alphatopo(:,s);
    ft_topoplotER(cfgtopo, topodata);
    
    %subplot(length(allsubj)+1,3,spos(sp)+1), hold on
    subplot(4,3,spos(sp)+1), hold on
    topodata.avg = betatopo(:,s);
    ft_topoplotER(cfgtopo, topodata);
    
    %subplot(length(allsubj)+1,3,spos(sp)+2), hold on
    subplot(4,3,spos(sp)+2), hold on
    topodata.avg = alltopo(:,s);
    ft_topoplotER(cfgtopo, topodata);
end

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- 2) PLOT CHOICE-PREDICTIVE LATERALIZATION OVER TIME/FREQ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and process data
TchoiceF_LI_all=[];
for s = 1:length(allsubj)
    load([loadpath,allsubj{s},'_trialwise_output_LIonly.mat'])
    TchoiceF_LI_all(:,:,s) = TchoiceF_LI(freqs>=1,:);
end
TchoiceF_LI_all(:,:,end+1) = mean(TchoiceF_LI_all,3);

% Plot
onsets = 0.4:0.4:(0.4*12);
snfig = [5 9 13 17];
sp=0;
figure, 
for s = 1:length(allsubj)+1
    if ~isempty(find(snfig==s))
        figure, sp=1;
    else
        sp=sp+1;
    end
    %subplot(2,length(allsubj)+1,sp), hold on
    subplot(2,4,sp), hold on
    imagesc(trltimes,freqs(freqs>=1),TchoiceF_LI_all(:,:,s),[-max(max(abs(TchoiceF_LI_all(:,:,s)))) max(max(abs(TchoiceF_LI_all(:,:,s))))].*0.65),
    line([0 0],[min(freqs(freqs>=1)) max(freqs)],'LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
    line([onsets; onsets],[repmat(min(freqs(freqs>=1)),1,length(onsets)); repmat(max(freqs(freqs>=1)),1,length(onsets))],'LineStyle','--','Color',[0.3 0.3 0.3])
    set(gca,'YDir','normal'), xlim([min(trltimes) max(trltimes)]); ylim([min(freqs(freqs>=1)) max(freqs)])
    if s==1, ylabel('Frequency (Hz)'), end
    if s==2, title('CHOICE-PREDICTIVE LATERALIZATION'), end    
    %subplot(2,length(allsubj)+1,sp+4), hold on
    subplot(2,4,sp+4), hold on
    plot(trltimes,mean(TchoiceF_LI_all(freqs>=alphabnd(1) & freqs<=alphabnd(2),:,s),1),'g')  % alpha = green
    plot(trltimes,mean(TchoiceF_LI_all(freqs>=betabnd(1) & freqs<=betabnd(2),:,s),1),'b')  % beta = blue
    plot(trltimes,mean(TchoiceF_LI_all(freqs>=allbnd(1) & freqs<=allbnd(2),:,s),1),'k')  % combined = black
    ylim([min(min(min(TchoiceF_LI_all))) max(max(max(TchoiceF_LI_all)))])
    line([0 0],get(gca, 'ylim')','LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
    line([onsets; onsets],repmat(get(gca, 'ylim')',1,length(onsets)),'LineStyle','--','Color',[0.5 0.5 0.5])
    line([min(trltimes) max(trltimes)],[0 0],'LineStyle','-','Color',[0.3 0.3 0.3])
    xlim([min(trltimes) max(trltimes)])
    if s==1, xlabel('Time relative to pre-mask (s)'), ylabel('T-score'), end
    if s<length(allsubj)+1; title(allsubj{s}); else title('MEAN'); end
end

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %-- 3) PLOT FINAL LPR-PREDICTIVE LATERALIZATION OVER TIME/FREQ
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and process data
TposteriorF_LI_all=[];
for s = 1:length(allsubj)
    load([loadpath,allsubj{s},'_trialwise_output_LIonly.mat'])
    TposteriorF_LI_all(:,:,s) = TposteriorF_LI(freqs>=1,:);
end
TposteriorF_LI_all(:,:,end+1) = mean(TposteriorF_LI_all,3);

% Plot
onsets = 0.4:0.4:(0.4*12);
snfig = [5 9 13 17];
sp=0;
figure,
for s = 1:length(allsubj)+1
    if ~isempty(find(snfig==s))
        figure, sp=1;
    else
        sp=sp+1;
    end
    subplot(2,4,sp), hold on
    imagesc(trltimes,freqs(freqs>=1),TposteriorF_LI_all(:,:,s),[-max(max(abs(TposteriorF_LI_all(:,:,s)))) max(max(abs(TposteriorF_LI_all(:,:,s))))].*0.65),
    line([0 0],[min(freqs(freqs>=1)) max(freqs)],'LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
    line([onsets; onsets],[repmat(min(freqs(freqs>=1)),1,length(onsets)); repmat(max(freqs(freqs>=1)),1,length(onsets))],'LineStyle','--','Color',[0.3 0.3 0.3])
    set(gca,'YDir','normal'), xlim([min(trltimes) max(trltimes)]); ylim([min(freqs(freqs>=1)) max(freqs)])
    if s==1, ylabel('Frequency (Hz)'), end
    if s==2, title('FINAL LPR-PREDICTIVE LATERALIZATION'), end
    
    subplot(2,4,sp+4), hold on
    plot(trltimes,mean(TposteriorF_LI_all(freqs>=alphabnd(1) & freqs<=alphabnd(2),:,s),1),'g')  % alpha = green
    plot(trltimes,mean(TposteriorF_LI_all(freqs>=betabnd(1) & freqs<=betabnd(2),:,s),1),'b')  % beta = blue
    plot(trltimes,mean(TposteriorF_LI_all(freqs>=allbnd(1) & freqs<=allbnd(2),:,s),1),'k')  % combined = black
    ylim([min(min(min(TposteriorF_LI_all))) max(max(max(TposteriorF_LI_all)))])
    line([0 0],get(gca, 'ylim')','LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
    line([onsets; onsets],repmat(get(gca, 'ylim')',1,length(onsets)),'LineStyle','--','Color',[0.5 0.5 0.5])
    line([min(trltimes) max(trltimes)],[0 0],'LineStyle','-','Color',[0.3 0.3 0.3])
    xlim([min(trltimes) max(trltimes)])
    if s==1, xlabel('Time relative to pre-mask (s)'), ylabel('T-score'), end
    if s<length(allsubj)+1; title(allsubj{s}); else title('MEAN'); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- 4) PLOT SAMPLE-ALIGNED LPR-PREDICTIVE LATERALIZATION OVER TIME/FREQ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_name = '_samplewise_output_LIonly.mat';
% Load and process data
TposteriorS_LI_all=[];
for s = 1:length(allsubj)
    load([loadpath,allsubj{s},file_name])
    TposteriorS_LI_all(:,:,:,s) = TposteriorS_LI(freqs>=1,:,:);
end
TposteriorS_LI_all(:,:,:,end+1) = nanmean(TposteriorS_LI_all,4);

% Plot
onsets = 0.4:0.4:(0.4*12);
gcols = [[0.6:-0.1:0 zeros(1,5)]; [ones(1,7) 0.9:-0.1:0.5]; [0.6:-0.1:0 zeros(1,5)]]';
bcols = [[0.6:-0.1:0 zeros(1,5)]; [0.6:-0.1:0 zeros(1,5)]; [ones(1,7) 0.9:-0.1:0.5]]';
kcols = [linspace(0.7,0,12); linspace(0.7,0,12); linspace(0.7,0,12)]';
snfig = [5 9 13 17];
sp=0;
figure,
for s = 1:length(allsubj)+1
    if ~isempty(find(snfig==s))
        figure, sp=1;
    else
        sp=sp+1;
    end
    subplot(2,4,sp), hold on
    imagesc(smptimes,freqs(freqs>=1),mean(TposteriorS_LI_all(:,:,:,s),3),[-max(max(abs(mean(TposteriorS_LI_all(:,:,:,s),3)))) max(max(abs(mean(TposteriorS_LI_all(:,:,:,s),3))))].*0.95),
    line([0 0],[min(freqs(freqs>=1)) max(freqs)],'LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
    set(gca,'YDir','normal'), xlim([min(smptimes) max(smptimes)]); ylim([min(freqs(freqs>=1)) max(freqs)])
    if s==1, xlabel('Time relative to sample onset (s)'), ylabel('Frequency (Hz)'), end
    if s==2, title('SAMPLE-ALIGNED LPR-PREDICTIVE LATERALIZATION'), end
        
    subplot(2,4,sp+4), hold on
    for samp = 1:length(onsets)
        plot(smptimes+onsets(samp),mean(TposteriorS_LI_all(freqs>=alphabnd(1) & freqs<=alphabnd(2),:,samp,s),1),'Color',gcols(samp,:))  % alpha = green
        plot(smptimes+onsets(samp),mean(TposteriorS_LI_all(freqs>=betabnd(1) & freqs<=betabnd(2),:,samp,s),1),'Color',bcols(samp,:))  % beta = blue
        plot(smptimes+onsets(samp),mean(TposteriorS_LI_all(freqs>=allbnd(1) & freqs<=allbnd(2),:,samp,s),1),'Color',kcols(samp,:))  % combined = black
    end
    ylim([min(min(min(min(TposteriorS_LI_all)))) max(max(max(max(TposteriorS_LI_all))))])
    line([0 0],get(gca, 'ylim')','LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
    line([onsets; onsets],repmat(get(gca, 'ylim')',1,length(onsets)),'LineStyle','--','Color',[0.5 0.5 0.5])
    line([-0.1 onsets(end)+smptimes(end)],[0 0],'LineStyle','-','Color',[0.3 0.3 0.3])
    xlim([-0.1 onsets(end)+smptimes(end)])
    if s==1, xlabel('Time relative to pre-mask (s)'), ylabel('T-score'), end
    if s<length(allsubj)+1; title(allsubj{s}); else title('MEAN'); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- 5) PLOT SAMPLE-ALIGNED LLR-PREDICTIVE LATERALIZATION OVER TIME/FREQ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and process data
file_name = '_samplewise_output_LIonly.mat';
TllrS_LI_all=[];
for s = 1:length(allsubj)
    load([loadpath,allsubj{s},file_name])
    TllrS_LI_all(:,:,:,s) = TllrS_LI(freqs>=1,:,:);
end
TllrS_LI_all(:,:,:,end+1) = nanmean(TllrS_LI_all,4);

% Plot
onsets = 0.4:0.4:(0.4*12);
gcols = [[0.6:-0.1:0 zeros(1,5)]; [ones(1,7) 0.9:-0.1:0.5]; [0.6:-0.1:0 zeros(1,5)]]';
bcols = [[0.6:-0.1:0 zeros(1,5)]; [0.6:-0.1:0 zeros(1,5)]; [ones(1,7) 0.9:-0.1:0.5]]';
kcols = [linspace(0.7,0,12); linspace(0.7,0,12); linspace(0.7,0,12)]';
snfig = [5 9 13 17];
sp=0;
figure, 
for s = 1:length(allsubj)+1
    if ~isempty(find(snfig==s))
        figure, sp=1;
    else
        sp=sp+1;
    end
    subplot(2,4,sp), hold on
    imagesc(smptimes,freqs(freqs>=1),mean(TllrS_LI_all(:,:,:,s),3),[-max(max(abs(mean(TllrS_LI_all(:,:,:,s),3)))) max(max(abs(mean(TllrS_LI_all(:,:,:,s),3))))].*0.95),
    line([0 0],[min(freqs(freqs>=1)) max(freqs)],'LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
    set(gca,'YDir','normal'), xlim([min(smptimes) max(smptimes)]); ylim([min(freqs(freqs>=1)) max(freqs)])
    if s==1, xlabel('Time relative to sample onset (s)'), ylabel('Frequency (Hz)'), end
    if s==2, title('SAMPLE-ALIGNED LLR-PREDICTIVE LATERALIZATION'), end
    
    subplot(2,4,sp+4), hold on
    for samp = 1:length(onsets)
        plot(smptimes+onsets(samp),mean(TllrS_LI_all(freqs>=alphabnd(1) & freqs<=alphabnd(2),:,samp,s),1),'Color',gcols(samp,:))  % alpha = green
        plot(smptimes+onsets(samp),mean(TllrS_LI_all(freqs>=betabnd(1) & freqs<=betabnd(2),:,samp,s),1),'Color',bcols(samp,:))  % beta = blue
        plot(smptimes+onsets(samp),mean(TllrS_LI_all(freqs>=allbnd(1) & freqs<=allbnd(2),:,samp,s),1),'Color',kcols(samp,:))  % combined = black
    end
    ylim([min(min(min(min(TllrS_LI_all)))) max(max(max(max(TllrS_LI_all))))])
    line([0 0],get(gca, 'ylim')','LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
    line([onsets; onsets],repmat(get(gca, 'ylim')',1,length(onsets)),'LineStyle','--','Color',[0.5 0.5 0.5])
    line([-0.1 onsets(end)+smptimes(end)],[0 0],'LineStyle','-','Color',[0.3 0.3 0.3])
    xlim([-0.1 onsets(end)+smptimes(end)])
    if s==1, xlabel('Time relative to pre-mask (s)'), ylabel('T-score'), end
    if s<length(allsubj)+1; title(allsubj{s}); else title('MEAN'); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- 6) PLOT SAMPLE-ALIGNED LLR*SURPRISE INTERACTION ON LATERALIZATION OVER TIME/FREQ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and process data
file_name = '_samplewise_output_LIonly.mat';
TpriorWllrS_LI_all=[]; TpriorllrS_LI_all=[]; TllrXsurpriseS_LI_all=[];
TprecllrS_LI_all=[];  TprecllrXsurpriseS_LI_all=[];  TprecWllrS_LI_all=[];
for s = 1:length(allsubj)
    load([loadpath,allsubj{s},file_name])
    TpriorWllrS_LI_all(:,:,:,s) = TpriorWllrS_LI(freqs>=1,:,:);
    TpriorllrS_LI_all(:,:,:,s) = TpriorllrS_LI(freqs>=1,:,:);
    TllrXsurpriseS_LI_all(:,:,:,s) = TllrXsurpriseS_LI(freqs>=1,:,:);
    
    TprecllrS_LI_all(:,:,:,s) = TprecllrS_LI(freqs>=1,:,:);
    TprecllrXsurpriseS_LI_all(:,:,:,s) = TprecllrXsurpriseS_LI(freqs>=1,:,:);
    TprecWllrS_LI_all(:,:,:,s) = TprecWllrS_LI(freqs>=1,:,:);   
end
TpriorWllrS_LI_all(:,:,:,end+1) = nanmean(TpriorWllrS_LI_all,4);
TpriorllrS_LI_all(:,:,:,end+1) = nanmean(TpriorllrS_LI_all,4);
TllrXsurpriseS_LI_all(:,:,:,end+1) = nanmean(TllrXsurpriseS_LI_all,4);

TprecllrS_LI_all(:,:,:,end+1) = nanmean(TprecllrS_LI_all,4);
TprecllrXsurpriseS_LI_all(:,:,:,end+1) = nanmean(TprecllrXsurpriseS_LI_all,4);
TprecWllrS_LI_all(:,:,:,end+1) = nanmean(TprecWllrS_LI_all,4);

% Plot
onsets = 0.4:0.4:(0.4*12);
gcols = [[0.6:-0.1:0 zeros(1,5)]; [ones(1,7) 0.9:-0.1:0.5]; [0.6:-0.1:0 zeros(1,5)]]';
bcols = [[0.6:-0.1:0 zeros(1,5)]; [0.6:-0.1:0 zeros(1,5)]; [ones(1,7) 0.9:-0.1:0.5]]';
kcols = [linspace(0.7,0,12); linspace(0.7,0,12); linspace(0.7,0,12)]';
snfig = [5 9 13 17];
sp=0;
figure,
for s = 1:length(allsubj)+1
    if ~isempty(find(snfig==s))
        figure, sp=1;
    else
        sp=sp+1;
    end
    subplot(4,4,sp), hold on
    imagesc(smptimes,freqs(freqs>=1),nanmean(TpriorWllrS_LI_all(:,:,:,s),3),[-max(max(abs(nanmean(TpriorWllrS_LI_all(:,:,:,s),3)))) max(max(abs(nanmean(TpriorWllrS_LI_all(:,:,:,s),3))))].*0.95),
    line([0 0],[min(freqs(freqs>=1)) max(freqs)],'LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
    set(gca,'YDir','normal'), xlim([min(smptimes) max(smptimes)]); ylim([min(freqs(freqs>=1)) max(freqs)])
    if s==1, xlabel('Time relative to sample onset (s)'), ylabel('Frequency (Hz)'), end
    if s==2, title('SAMPLE-ALIGNED LLR*SURPRISE INTERACTION ON LATERALIZATION'), end
    
    subplot(4,4,sp+4), hold on
    imagesc(smptimes,freqs(freqs>=1),nanmean(TpriorllrS_LI_all(:,:,:,s),3),[-max(max(abs(nanmean(TpriorllrS_LI_all(:,:,:,s),3)))) max(max(abs(nanmean(TpriorllrS_LI_all(:,:,:,s),3))))].*0.95),
    line([0 0],[min(freqs(freqs>=1)) max(freqs)],'LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
    set(gca,'YDir','normal'), xlim([min(smptimes) max(smptimes)]); ylim([min(freqs(freqs>=1)) max(freqs)])
    if s==1, xlabel('Time relative to sample onset (s)'), ylabel('Frequency (Hz)'), end
    
    subplot(4,4,sp+8), hold on
    imagesc(smptimes,freqs(freqs>=1),nanmean(TllrXsurpriseS_LI_all(:,:,:,s),3),[-max(max(abs(nanmean(TllrXsurpriseS_LI_all(:,:,:,s),3)))) max(max(abs(nanmean(TllrXsurpriseS_LI_all(:,:,:,s),3))))].*0.95),
    line([0 0],[min(freqs(freqs>=1)) max(freqs)],'LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
    set(gca,'YDir','normal'), xlim([min(smptimes) max(smptimes)]); ylim([min(freqs(freqs>=1)) max(freqs)])
    if s==1, xlabel('Time relative to sample onset (s)'), ylabel('Frequency (Hz)'), end
    
    subplot(4,4,sp+12), hold on
    for samp = 2:length(onsets)
        plot(smptimes+onsets(samp),mean(TllrXsurpriseS_LI_all(freqs>=alphabnd(1) & freqs<=alphabnd(2),:,samp-1,s),1),'Color',gcols(samp,:))  % alpha = green
        plot(smptimes+onsets(samp),mean(TllrXsurpriseS_LI_all(freqs>=betabnd(1) & freqs<=betabnd(2),:,samp-1,s),1),'Color',bcols(samp,:))  % beta = blue
        plot(smptimes+onsets(samp),mean(TllrXsurpriseS_LI_all(freqs>=allbnd(1) & freqs<=allbnd(2),:,samp-1,s),1),'Color',kcols(samp,:))  % combined = black
    end
    if s<length(allsubj)+1, ylim([min(min(min(min(TllrXsurpriseS_LI_all)))) max(max(max(max(TllrXsurpriseS_LI_all))))]), end
    line([0 0],get(gca, 'ylim')','LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
    line([onsets; onsets],repmat(get(gca, 'ylim')',1,length(onsets)),'LineStyle','--','Color',[0.5 0.5 0.5])
    line([-0.1 onsets(end)+smptimes(end)],[0 0],'LineStyle','-','Color',[0.3 0.3 0.3])
    xlim([-0.1 onsets(end)+smptimes(end)])
    if s==1, xlabel('Time relative to pre-mask (s)'), ylabel('T-score'), end
    if s<length(allsubj)+1; title(allsubj{s}); else title('MEAN'); end
    if s<length(allsubj)+1; title(allsubj{s}); else title('MEAN'); end
end

% Plot precedent timestamp power approach
sp=0;
figure,
for s = 1:length(allsubj)+1
    if ~isempty(find(snfig==s))
        figure, sp=1;
    else
        sp=sp+1;
    end
    subplot(4,4,sp), hold on
 
    imagesc(smptimes,freqs(freqs>=1),nanmean(TprecWllrS_LI_all(:,:,:,s),3),[-max(max(abs(nanmean(TprecWllrS_LI_all(:,:,:,s),3)))) max(max(abs(nanmean(TprecWllrS_LI_all(:,:,:,s),3))))].*0.95),
    line([0 0],[min(freqs(freqs>=1)) max(freqs)],'LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
    set(gca,'YDir','normal'), xlim([min(smptimes) max(smptimes)]); ylim([min(freqs(freqs>=1)) max(freqs)])
    if s==1, xlabel('Time relative to sample onset (s)'), ylabel('Frequency (Hz)'), end
    if s==2, title('SAMPLE-ALIGNED LLR*SURPRISE INTERACTION ON LATERALIZATION (Precedent timestamp power)'), end
    
    subplot(4,4,sp+4), hold on
    imagesc(smptimes,freqs(freqs>=1),nanmean(TprecllrS_LI_all(:,:,:,s),3),[-max(max(abs(nanmean(TprecllrS_LI_all(:,:,:,s),3)))) max(max(abs(nanmean(TprecllrS_LI_all(:,:,:,s),3))))].*0.95),
    line([0 0],[min(freqs(freqs>=1)) max(freqs)],'LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
    set(gca,'YDir','normal'), xlim([min(smptimes) max(smptimes)]); ylim([min(freqs(freqs>=1)) max(freqs)])
    if s==1, xlabel('Time relative to sample onset (s)'), ylabel('Frequency (Hz)'), end
    
    subplot(4,4,sp+8), hold on
    imagesc(smptimes,freqs(freqs>=1),nanmean(TprecllrXsurpriseS_LI_all(:,:,:,s),3),[-max(max(abs(nanmean(TprecllrXsurpriseS_LI_all(:,:,:,s),3)))) max(max(abs(nanmean(TprecllrXsurpriseS_LI_all(:,:,:,s),3))))].*0.95),
    line([0 0],[min(freqs(freqs>=1)) max(freqs)],'LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
    set(gca,'YDir','normal'), xlim([min(smptimes) max(smptimes)]); ylim([min(freqs(freqs>=1)) max(freqs)])
    if s==1, xlabel('Time relative to sample onset (s)'), ylabel('Frequency (Hz)'), end
    
    subplot(4,4,sp+12), hold on
    for samp = 2:length(onsets)
        plot(smptimes+onsets(samp),mean(TprecllrXsurpriseS_LI_all(freqs>=alphabnd(1) & freqs<=alphabnd(2),:,samp-1,s),1),'Color',gcols(samp,:))  % alpha = green
        plot(smptimes+onsets(samp),mean(TprecllrXsurpriseS_LI_all(freqs>=betabnd(1) & freqs<=betabnd(2),:,samp-1,s),1),'Color',bcols(samp,:))  % beta = blue
        plot(smptimes+onsets(samp),mean(TprecllrXsurpriseS_LI_all(freqs>=allbnd(1) & freqs<=allbnd(2),:,samp-1,s),1),'Color',kcols(samp,:))  % combined = black
    end
    if s<length(allsubj)+1, ylim([min(min(min(min(TprecllrXsurpriseS_LI_all)))) max(max(max(max(TprecllrXsurpriseS_LI_all))))]), end
    line([0 0],get(gca, 'ylim')','LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
    line([onsets; onsets],repmat(get(gca, 'ylim')',1,length(onsets)),'LineStyle','--','Color',[0.5 0.5 0.5])
    line([-0.1 onsets(end)+smptimes(end)],[0 0],'LineStyle','-','Color',[0.3 0.3 0.3])
    xlim([-0.1 onsets(end)+smptimes(end)])
    if s==1, xlabel('Time relative to pre-mask (s)'), ylabel('T-score'), end
    if s<length(allsubj)+1; title(allsubj{s}); else title('MEAN'); end
    if s<length(allsubj)+1; title(allsubj{s}); else title('MEAN'); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- 7) PLOT TRIAL-BY-TRIAL RELATIONSHIP B/W FINAL PRIOR & LATE LATERALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nbins = 20;
% Load and process data
fLI_alpha=[]; fLI_beta=[]; fLI_all=[]; f_prior=[];
figure,
for s = 1:length(allsubj)
    load([loadpath,allsubj{s},'_trialwise_output_LIonly.mat'])
    
    subplot(3,length(allsubj)+1,s), hold on, plot_binned_scatter(final_prior,mean(finalLI(:,freqs>=alphabnd(1) & freqs<=alphabnd(2)),2),nbins,'none',[0 1 0])
    r = corr(final_prior,mean(finalLI(:,freqs>=alphabnd(1) & freqs<=alphabnd(2)),2),'type','Spearman'); title(['rho = ',num2str(round(r,3))]), set(gca,'TickDir','out')
    subplot(3,length(allsubj)+1,s+4), hold on, plot_binned_scatter(final_prior,mean(finalLI(:,freqs>=betabnd(1) & freqs<=betabnd(2)),2),nbins,'none',[0 0 1])
    r = corr(final_prior,mean(finalLI(:,freqs>=betabnd(1) & freqs<=betabnd(2)),2),'type','Spearman'); title(['rho = ',num2str(round(r,3))]), set(gca,'TickDir','out')
    subplot(3,length(allsubj)+1,s+8), hold on, plot_binned_scatter(final_prior,mean(finalLI(:,freqs>=allbnd(1) & freqs<=allbnd(2)),2),nbins,'none',[0 0 0])
    r = corr(final_prior,mean(finalLI(:,freqs>=allbnd(1) & freqs<=allbnd(2)),2),'type','Spearman'); title(['rho = ',num2str(round(r,3))]), set(gca,'TickDir','out')
    if s==1, xlabel('Prior belief'), ylabel('Laterlization index'), end
    if s==2, title('TRIAL-BY-TRIAL RELATIONSHIP B/W FINAL PRIOR & LATE LATERALIZATION'), end
    
    fLI_alpha = [fLI_alpha; mean(finalLI(:,freqs>=alphabnd(1) & freqs<=alphabnd(2)),2)];
    fLI_beta = [fLI_beta; mean(finalLI(:,freqs>=betabnd(1) & freqs<=betabnd(2)),2)];
    fLI_all = [fLI_all; mean(finalLI(:,freqs>=allbnd(1) & freqs<=allbnd(2)),2)];
    f_prior = [f_prior; final_prior];
end
subplot(3,length(allsubj)+1,s+1), hold on, plot_binned_scatter(f_prior,fLI_alpha,nbins,'none',[0 1 0])
r = corr(f_prior,fLI_alpha,'type','Spearman'); title(['rho = ',num2str(round(r,3))]), set(gca,'TickDir','out')
subplot(3,length(allsubj)+1,s+1+4), hold on, plot_binned_scatter(f_prior,fLI_beta,nbins,'none',[0 0 1])
r = corr(f_prior,fLI_beta,'type','Spearman'); title(['rho = ',num2str(round(r,3))]), set(gca,'TickDir','out')
subplot(3,length(allsubj)+1,s+1+8), hold on, plot_binned_scatter(f_prior,fLI_all,nbins,'none',[0 0 0])
r = corr(f_prior,fLI_all,'type','Spearman'); title(['rho = ',num2str(round(r,3))]), set(gca,'TickDir','out')
   
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   STACKED SAMPLES APPROACH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_name = '_samplewise_stacksam_LIonly.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- 8) PLOT SAMPLE-ALIGNED LPR-PREDICTIVE LATERALIZATION OVER TIME/FREQ - STACKED APPROACH
%       (Analogous to 4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_name = '_samplewise_stacksam_LIonly.mat';
% Load and process data
TposteriorS_LI_all=[];
for s = 1:length(allsubj)
    load([loadpath,allsubj{s},file_name])
    TposteriorS_LI_all(:,:,s) = TposteriorS_LI(freqs>=1,:);
end
TposteriorS_LI_all(:,:,end+1) = mean(TposteriorS_LI_all,3);

% Plot
% onsets = 0:0.2:0.8;%0.4:0.4:(0.4*12);
gcols = [[0.6:-0.1:0 zeros(1,5)]; [ones(1,7) 0.9:-0.1:0.5]; [0.6:-0.1:0 zeros(1,5)]]';
bcols = [[0.6:-0.1:0 zeros(1,5)]; [0.6:-0.1:0 zeros(1,5)]; [ones(1,7) 0.9:-0.1:0.5]]';
kcols = [linspace(0.7,0,12); linspace(0.7,0,12); linspace(0.7,0,12)]';
snfig = [5 9 13 17];
sp=0;
figure,
for s = 1:length(allsubj)+1
    if ~isempty(find(snfig==s))
        figure, sp=1;
    else
        sp=sp+1;
    end
    subplot(2,4,sp), hold on
    imagesc(smptimes,freqs(freqs>=1),TposteriorS_LI_all(:,:,s),[-max(max(abs(TposteriorS_LI_all(:,:,s)))) max(max(abs(TposteriorS_LI_all(:,:,s))))].*0.65),
    line([0 0],[min(freqs(freqs>=1)) max(freqs)],'LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
    set(gca,'YDir','normal'), xlim([min(smptimes) max(smptimes)]); ylim([min(freqs(freqs>=1)) max(freqs)])
    if s==1, xlabel('Time relative to sample onset (s)'), ylabel('Frequency (Hz)'), end
    if s==2, title('SAMPLE-ALIGNED LPR-PREDICTIVE LATERALIZATION - STACKED APPROACH'), end
    
    subplot(2,4,sp+4), hold on
    plot(smptimes,mean(TposteriorS_LI_all(freqs>=alphabnd(1) & freqs<=alphabnd(2),:,s),1),'Color',gcols(2,:))  % alpha = green
    plot(smptimes,mean(TposteriorS_LI_all(freqs>=betabnd(1) & freqs<=betabnd(2),:,s),1),'Color',bcols(2,:))  % beta = blue
    plot(smptimes,mean(TposteriorS_LI_all(freqs>=allbnd(1) & freqs<=allbnd(2),:,s),1),'Color',kcols(2,:))  % combined = black

    ylim([min(min(min(TposteriorS_LI_all))) max(max(max(TposteriorS_LI_all)))])
    line([0 0],get(gca, 'ylim')','LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
%     line([onsets; onsets],repmat(get(gca, 'ylim')',1,length(onsets)),'LineStyle','--','Color',[0.5 0.5 0.5])
%     line([-0.1 onsets(end)],[0 0],'LineStyle','-','Color',[0.3 0.3 0.3])
%     xlim([-0.1 0.8])
    if s==1, xlabel('Time relative to pre-mask (s)'), ylabel('T-score'), end
    if s<length(allsubj)+1; title(allsubj{s}); else title('MEAN'); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- 9) PLOT SAMPLE-ALIGNED LLR-PREDICTIVE LATERALIZATION OVER TIME/FREQ -STACKED APPROACH
%--     (analogous to: 5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load and process data
file_name = '_samplewise_stacksam_LIonly.mat';
TllrS_LI_all=[];
for s = 1:length(allsubj)
    load([loadpath,allsubj{s},file_name])
    TllrS_LI_all(:,:,s) = TllrS_LI(freqs>=1,:);
end
TllrS_LI_all(:,:,end+1) = mean(TllrS_LI_all,3);

% Plot
% onsets = 0:0.2:0.8;%0.4:0.4:(0.4*12);
gcols = [[0.6:-0.1:0 zeros(1,5)]; [ones(1,7) 0.9:-0.1:0.5]; [0.6:-0.1:0 zeros(1,5)]]';
bcols = [[0.6:-0.1:0 zeros(1,5)]; [0.6:-0.1:0 zeros(1,5)]; [ones(1,7) 0.9:-0.1:0.5]]';
kcols = [linspace(0.7,0,12); linspace(0.7,0,12); linspace(0.7,0,12)]';
snfig = [5 9 13 17];
sp=0;
figure, 
for s = 1:length(allsubj)+1
    if ~isempty(find(snfig==s))
        figure, sp=1;
    else
        sp=sp+1;
    end
    subplot(2,4,sp), hold on
    imagesc(smptimes,freqs(freqs>=1),TllrS_LI_all(:,:,s),[-max(max(abs(TllrS_LI_all(:,:,s)))) max(max(abs(TllrS_LI_all(:,:,s))))].*0.65),
    line([0 0],[min(freqs(freqs>=1)) max(freqs)],'LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
    set(gca,'YDir','normal'), xlim([min(smptimes) max(smptimes)]); ylim([min(freqs(freqs>=1)) max(freqs)])
    if s==1, xlabel('Time relative to sample onset (s)'), ylabel('Frequency (Hz)'), end
    if s==2, title('SAMPLE-ALIGNED LLR-PREDICTIVE LATERALIZATION -STACKED APPROACH'), end
    
    subplot(2,4,sp+4), hold on
    plot(smptimes,mean(TllrS_LI_all(freqs>=alphabnd(1) & freqs<=alphabnd(2),:,s),1),'Color',gcols(2,:))  % alpha = green
    plot(smptimes,mean(TllrS_LI_all(freqs>=betabnd(1) & freqs<=betabnd(2),:,s),1),'Color',bcols(2,:))  % beta = blue
    plot(smptimes,mean(TllrS_LI_all(freqs>=allbnd(1) & freqs<=allbnd(2),:,s),1),'Color',kcols(2,:))  % combined = black

    ylim([min(min(min(TllrS_LI_all))) max(max(max(TllrS_LI_all)))])
    line([0 0],get(gca, 'ylim')','LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
%     line([onsets; onsets],repmat(get(gca, 'ylim')',1,length(onsets)),'LineStyle','--','Color',[0.5 0.5 0.5])
%     line([-0.1 onsets(end)],[0 0],'LineStyle','-','Color',[0.3 0.3 0.3])
%     xlim([-0.1 0.8])
    if s==1, xlabel('Time relative to pre-mask (s)'), ylabel('T-score'), end
    if s<length(allsubj)+1; title(allsubj{s}); else title('MEAN'); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- 10) PLOT SAMPLE-ALIGNED LLR*SURPRISE INTERACTION ON LATERALIZATION OVER TIME/FREQ
%   STACKED APPROACH
%   (ANALOGOUS TO 6)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load and process data
file_name = '_samplewise_stacksam_LIonly.mat';
TpriorWllrS_LI_all=[]; TpriorllrS_LI_all=[]; TllrXsurpriseS_LI_all=[]; 
TprecllrS_LI_all=[]; TprecllrXsurpriseS_LI_all=[]; TprecWllrS_LI_all=[];
for s = 1:length(allsubj)
    load([loadpath,allsubj{s},file_name])
    TpriorWllrS_LI_all(:,:,s) = TpriorWllrS_LI(freqs>=1,:);
    TpriorllrS_LI_all(:,:,s) = TpriorllrS_LI(freqs>=1,:);
    TllrXsurpriseS_LI_all(:,:,s) = TllrXsurpriseS_LI(freqs>=1,:);
    TprecllrS_LI_all(:,:,s) = TprecllrS_LI(freqs>=1,:);
    TprecllrXsurpriseS_LI_all(:,:,s) = TprecllrXsurpriseS_LI(freqs>=1,:);
    TprecWllrS_LI_all(:,:,s) = TprecWllrS_LI(freqs>=1,:);   
end
TpriorWllrS_LI_all(:,:,end+1) = nanmean(TpriorWllrS_LI_all,3);
TpriorllrS_LI_all(:,:,end+1) = nanmean(TpriorllrS_LI_all,3);
TllrXsurpriseS_LI_all(:,:,end+1) = nanmean(TllrXsurpriseS_LI_all,3);

TprecllrS_LI_all(:,:,end+1) = nanmean(TprecllrS_LI_all,3);
TprecllrXsurpriseS_LI_all(:,:,end+1) = nanmean(TprecllrXsurpriseS_LI_all,3);
TprecWllrS_LI_all(:,:,end+1) = nanmean(TprecWllrS_LI_all,3);

% Plot
onsets = 0.4:0.4:(0.4*12);
gcols = [[0.6:-0.1:0 zeros(1,5)]; [ones(1,7) 0.9:-0.1:0.5]; [0.6:-0.1:0 zeros(1,5)]]';
bcols = [[0.6:-0.1:0 zeros(1,5)]; [0.6:-0.1:0 zeros(1,5)]; [ones(1,7) 0.9:-0.1:0.5]]';
kcols = [linspace(0.7,0,12); linspace(0.7,0,12); linspace(0.7,0,12)]';
snfig = [5 9 13 17];
sp=0;
figure,
for s = 1:length(allsubj)+1
    if ~isempty(find(snfig==s))
        figure, sp=1;
    else
        sp=sp+1;
    end
    subplot(4,4,sp), hold on
    imagesc(smptimes,freqs(freqs>=1),TpriorWllrS_LI_all(:,:,s),[-max(max(abs(TpriorWllrS_LI_all(:,:,s)))) max(max(abs(TpriorWllrS_LI_all(:,:,s))))].*0.95),
    line([0 0],[min(freqs(freqs>=1)) max(freqs)],'LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
    set(gca,'YDir','normal'), xlim([min(smptimes) max(smptimes)]); ylim([min(freqs(freqs>=1)) max(freqs)])
    if s==1, xlabel('Time relative to sample onset (s)'), ylabel('Frequency (Hz)'), end
    if s==2, title('SAMPLE-ALIGNED LLR*SURPRISE INTERACTION ON LATERALIZATION - STACKED APPROACH'), end
    
    subplot(4,4,sp+4), hold on
    imagesc(smptimes,freqs(freqs>=1),TpriorllrS_LI_all(:,:,s),[-max(max(abs(TpriorllrS_LI_all(:,:,s)))) max(max(abs(TpriorllrS_LI_all(:,:,s))))].*0.95),
    line([0 0],[min(freqs(freqs>=1)) max(freqs)],'LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
    set(gca,'YDir','normal'), xlim([min(smptimes) max(smptimes)]); ylim([min(freqs(freqs>=1)) max(freqs)])
    if s==1, xlabel('Time relative to sample onset (s)'), ylabel('Frequency (Hz)'), end
    
    subplot(4,4,sp+8), hold on
    imagesc(smptimes,freqs(freqs>=1),TllrXsurpriseS_LI_all(:,:,s),[-max(max(abs(TllrXsurpriseS_LI_all(:,:,s)))) max(max(abs(TllrXsurpriseS_LI_all(:,:,s))))].*0.95),
    line([0 0],[min(freqs(freqs>=1)) max(freqs)],'LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
    set(gca,'YDir','normal'), xlim([min(smptimes) max(smptimes)]); ylim([min(freqs(freqs>=1)) max(freqs)])
    if s==1, xlabel('Time relative to sample onset (s)'), ylabel('Frequency (Hz)'), end
    
    subplot(4,4,sp+12), hold on
    for samp = 2:length(onsets)
        plot(smptimes,mean(TllrXsurpriseS_LI_all(freqs>=alphabnd(1) & freqs<=alphabnd(2),:,s),1),'Color',gcols(samp,:))  % alpha = green
        plot(smptimes,mean(TllrXsurpriseS_LI_all(freqs>=betabnd(1) & freqs<=betabnd(2),:,s),1),'Color',bcols(samp,:))  % beta = blue
        plot(smptimes,mean(TllrXsurpriseS_LI_all(freqs>=allbnd(1) & freqs<=allbnd(2),:,s),1),'Color',kcols(samp,:))  % combined = black
    end
    if s<length(allsubj)+1, ylim([min(min(min(TllrXsurpriseS_LI_all))) max(max(max(TllrXsurpriseS_LI_all)))]), end
    line([0 0],get(gca, 'ylim')','LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
    line([onsets; onsets],repmat(get(gca, 'ylim')',1,length(onsets)),'LineStyle','--','Color',[0.5 0.5 0.5])
    line([-0.1 onsets(end)+smptimes(end)],[0 0],'LineStyle','-','Color',[0.3 0.3 0.3])
   xlim([-0.1 smptimes(end)])
    if s==1, xlabel('Time relative to pre-mask (s)'), ylabel('T-score'), end
    if s<length(allsubj)+1; title(allsubj{s}); else title('MEAN'); end
end

% Plot precedent time stamp power
sp=0;
figure,
for s = 1:length(allsubj)+1
    if ~isempty(find(snfig==s))
        figure, sp=1;
    else
        sp=sp+1;
    end
    subplot(4,4,sp), hold on
    imagesc(smptimes,freqs(freqs>=1),TprecWllrS_LI_all(:,:,s),[-max(max(abs(TprecWllrS_LI_all(:,:,s)))) max(max(abs(TprecWllrS_LI_all(:,:,s))))].*0.95),
    line([0 0],[min(freqs(freqs>=1)) max(freqs)],'LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
    set(gca,'YDir','normal'), xlim([min(smptimes) max(smptimes)]); ylim([min(freqs(freqs>=1)) max(freqs)])
    if s==1, xlabel('Time relative to sample onset (s)'), ylabel('Frequency (Hz)'), end
    if s==2, title('SAMPLE-ALIGNED LLR*SURPRISE INTERACTION ON LATERALIZATION - STACKED APPROACH'), end
    
    subplot(4,4,sp+4), hold on
    imagesc(smptimes,freqs(freqs>=1),TprecllrS_LI_all(:,:,s),[-max(max(abs(TprecllrS_LI_all(:,:,s)))) max(max(abs(TprecllrS_LI_all(:,:,s))))].*0.95),
    line([0 0],[min(freqs(freqs>=1)) max(freqs)],'LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
    set(gca,'YDir','normal'), xlim([min(smptimes) max(smptimes)]); ylim([min(freqs(freqs>=1)) max(freqs)])
    if s==1, xlabel('Time relative to sample onset (s)'), ylabel('Frequency (Hz)'), end
    
    subplot(4,4,sp+8), hold on
    imagesc(smptimes,freqs(freqs>=1),TprecllrXsurpriseS_LI_all(:,:,s),[-max(max(abs(TprecllrXsurpriseS_LI_all(:,:,s)))) max(max(abs(TprecllrXsurpriseS_LI_all(:,:,s))))].*0.95),
    line([0 0],[min(freqs(freqs>=1)) max(freqs)],'LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
    set(gca,'YDir','normal'), xlim([min(smptimes) max(smptimes)]); ylim([min(freqs(freqs>=1)) max(freqs)])
    if s==1, xlabel('Time relative to sample onset (s)'), ylabel('Frequency (Hz)'), end
    
    subplot(4,4,sp+12), hold on
    for samp = 2:length(onsets)
        plot(smptimes,mean(TprecllrXsurpriseS_LI_all(freqs>=alphabnd(1) & freqs<=alphabnd(2),:,s),1),'Color',gcols(samp,:))  % alpha = green
        plot(smptimes,mean(TprecllrXsurpriseS_LI_all(freqs>=betabnd(1) & freqs<=betabnd(2),:,s),1),'Color',bcols(samp,:))  % beta = blue
        plot(smptimes,mean(TprecllrXsurpriseS_LI_all(freqs>=allbnd(1) & freqs<=allbnd(2),:,s),1),'Color',kcols(samp,:))  % combined = black
    end
    if s<length(allsubj)+1, ylim([min(min(min(TprecllrXsurpriseS_LI_all))) max(max(max(TprecllrXsurpriseS_LI_all)))]), end
    line([0 0],get(gca, 'ylim')','LineStyle','-','Color',[0 0 0],'LineWidth',1.25)
    line([onsets; onsets],repmat(get(gca, 'ylim')',1,length(onsets)),'LineStyle','--','Color',[0.5 0.5 0.5])
    line([-0.1 onsets(end)+smptimes(end)],[0 0],'LineStyle','-','Color',[0.3 0.3 0.3])
   xlim([-0.1 smptimes(end)])
    if s==1, xlabel('Time relative to pre-mask (s)'), ylabel('T-score'), end
    if s<length(allsubj)+1; title(allsubj{s}); else title('MEAN'); end
end