clear

addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
addpath '/mnt/homes/home024/chernandez/matlab/cbrewer'
ft_defaults
loadpath = '/mnt/homes/home024/chernandez/meg_data/surprise/preprocessed/Data/ERF/';  % path of preprocessed MEG data

% load some stuff for plotting
load([loadpath,'CTF275_helmet.mat']);  % load layout
lay.outline = lay.outline([1 3:end]); % remove outer bound

[T,cmap] = evalc(['cbrewer(''div'', ''RdBu'', 64)']);  % create colourmap (and suppressing warning messages using evalc) - Caro, you need to download the cbrewer set of functions for this to work...
cmap = flipud(cmap); % red = increase                  % ... get those here: https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab

% model variables
allvars = {'TdeltaL2','Tposterior1','Tposterior2','Tprior1', 'Tprior2','Tsurprise2'};
%allvars = {'Tsurprise2'};
all_ch_vars = {'Tchoice1','Tchoice2'};

%subjects
allsubj = {'DHB','EXF','TFD'};
%allsubj = {'DHB'};

%Define the time windows
dot_times = 0:0.01:0.39;
choice_times = 0:0.01:0.59;
t_window = [0.05 0.150 0.225 0.3 0.4];
t_window_choice = [0.5 0.6];


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
cnt=1;

for i = 1:length(allsubj),%%Subjects
    subject = allsubj{i};
    load([loadpath,subject,'_regression_output.mat']);
    Tchoice1All(:,:,i) = Tchoice1;
    Tchoice2All(:,:,i) = Tchoice2;
    TdeltaL2All(:,:,:,i) = TdeltaL2;
    Tposterior1All(:,:,:,i) = Tposterior1;
    Tposterior2All(:,:,:,i) = Tposterior2;
    Tprior1All(:,:,:,i) = Tprior1;
    Tprior2All(:,:,:,i) = Tprior2;
    Tsurprise2All(:,:,:,i) = Tsurprise2;
   
end
%Mean all subjects
Tchoice1All(:,:,i+1) = mean(Tchoice1All(:,:,:),3);
Tchoice2All(:,:,i+1) = mean(Tchoice2All(:,:,:),3);
TdeltaL2All(:,:,:,i+1) = mean(TdeltaL2All(:,:,:,:),4);
Tposterior1All(:,:,:,i+1) = mean(Tposterior1All(:,:,:,:),4);
Tposterior2All(:,:,:,i+1) = mean(Tposterior2All(:,:,:,:),4);
Tprior1All(:,:,:,i+1) = mean(Tprior1All(:,:,:,:),4);
Tprior2All(:,:,:,i+1) = mean(Tprior2All(:,:,:,:),4);
Tsurprise2All(:,:,:,i+1) = mean(Tsurprise2All(:,:,:,:),4);

%Get rid of useless channels
%Here, assuming all grad are the same
grad.chanpos = grad.chanpos(1:274,:);
grad.chanori = grad.chanori(1:274,:);
grad.chantype = grad.chantype(1:274);
grad.chanunit = grad.chanunit(1:274);
grad.label = grad.label(1:274);

% Create data structure for plotting
topodata.dimord     = 'chan_time';
topodata.label      = grad.label;
topodata.time       = 0;
topodata.grad       = grad;

%Choice variables
for i=1:length(all_ch_vars), %%Model Variables
    cnt=1;
    figure,
    suptitle(all_ch_vars{i});
    for j = 1:length(allsubj)+1,%%Subjects and mean of subjects
        ini_w = 0.4;
        for k=1: length(t_window_choice),%Time windows
            times_w = (ini_w<=choice_times & choice_times<t_window_choice(k));
            ini_w=t_window_choice(k);
            times_w = find(times_w); %indexes of times of interest
            switch all_ch_vars{i}
                case 'Tchoice1'
                    topodata.avg = mean(mean(Tchoice1All(:,times_w,j),2),3);
                 case 'Tchoice2'
                     topodata.avg = mean(mean(Tchoice2All(:,times_w,j),2),3);
            end

            % Set colour range
            cpeak = max(abs(topodata.avg));
            cfgtopo.zlim        = [-cpeak(1) cpeak(1)];

            %TODO: fix subplot titles, these are not working
            s(cnt)=subplot(4,2,cnt);  
            if cnt<6  
                title(s(cnt),[ini_w '-' t_window_choice(k)]);
            end
            hold on,
            ft_topoplotER(cfgtopo, topodata);
            cnt = cnt + 1;
        end
    end
end            
for i=1:length(allvars), %%Model Variables
    cnt=1;
    figure,
    suptitle(allvars{i});
    for j = 1:length(allsubj)+1,%%Subjects and mean of subjects
        ini_w = 0;
        for k=1: length(t_window),%Time windows
            times_w = (ini_w<=dot_times & dot_times<t_window(k));
            ini_w=t_window(k);
            times_w = find(times_w); %indexes of times of interest
            switch allvars{i}
                case 'TdeltaL2'
                    topodata.avg = mean(mean(TdeltaL2All(:,:,times_w,j),3),2); %all dots
%                     topodata.avg = mean(TdeltaL2All(:,11,times_w,j),3); %last dot
                case 'Tposterior1'
                    topodata.avg= mean(mean(Tposterior1All(:,:,times_w,j),3),2); %all dots
%                     topodata.avg= mean(Tposterior1All(:,12,times_w,j),3),2); %last dot
                case 'Tposterior2'
                    topodata.avg= mean(mean(Tposterior2All(:,:,times_w,j),3),2); %all dots
%                     topodata.avg= mean(Tposterior2All(:,12,times_w,j),3),2); %last dot
                case 'Tprior1'
                    topodata.avg= mean(mean(Tprior1All(:,:,times_w,j),3),2); %all dots
%                     topodata.avg= mean(Tprior1All(:,11,times_w,j),3),2); %last dot
                case 'Tprior2'
                    topodata.avg= mean(mean(Tprior2All(:,:,times_w,j),3),2); %all dots
%                     topodata.avg= mean(Tprior2All(:,11,times_w,j),3),2); %last dot
                case 'Tsurprise2'
                    topodata.avg= mean(mean(Tsurprise2All(:,:,times_w,j),3),2); %all dots
%                     topodata.avg= mean(Tsurprise2All(:,11,times_w,j),3),2); %last dot
            end
            
            % Set colour range
            cpeak = max(abs(topodata.avg));
            cfgtopo.zlim        = [-cpeak(1) cpeak(1)];
            
            %TODO: fix subplot titles, these are not working
            s(cnt)=subplot(4,5,cnt);  
            if cnt<6  
                title(s(cnt),[ini_w '-' t_window(k)]);
            end
            hold on,
%               if cnt<6             
%                 title(s(cnt),[ini_w '-' t_window(k)]),
%                 grid on ,
%                 xlabel( 'window' ) ;  ylabel( 'subject' ) ,
%            end

            %title([allsubj{i},': ' ,timew(k,1),'-',timew(k,2),' ms']),
            ft_topoplotER(cfgtopo, topodata);
            cnt = cnt + 1;
        end
         
    end
end
 
