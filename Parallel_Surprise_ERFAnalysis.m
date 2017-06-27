% Top-level script for running multiple 2-accumulator LCA fits in parallel
% Uses qsubcellfun.m to submit jobs to TORQUE


% Path stuff
addpath(genpath('/mnt/homes/home024/chernandez'))
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'
ft_defaults
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221/qsub'
%addpath '/mnt/homes/home024/pmurphy/Surprise_accumulation/Analysis/MEG/Scripts/'

% Subject/model variant stuff
allsubj = {'EXF' 'TFD' 'DHB'};

% Submit jobs to TORQUE
setenv('TORQUEHOME', 'yes')   % not sure what this does..........   
cd('~/Data/')

memreq = 4*1024^3;  % memory required (1GB = 1024^3)
timreq = 1*60*60;   % time required (hours*minutes*seconds)

qsubcellfun(@runERFanalysis, allsubj, 'compile', 'no',...
         'memreq', memreq, 'timreq', timreq, 'stack', 1, 'StopOnError', false, 'backend', 'torque');
    
