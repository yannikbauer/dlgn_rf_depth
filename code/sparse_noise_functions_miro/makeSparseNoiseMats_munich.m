%% makeSparseNoiseMats_munich
% Makes RF maps from MUA-filtered recordings to the sparse noise stimulus
% NOTE: previous testSparseNoiseScript in obsolete and test
% TODO: 

%% Setup
clear all; close all;
startup_lmu;
% cd /Volumes/lab/users/miro/code/sparseNoiseFunctions
cd(fileparts(which(mfilename))); % changes working dir to current file
addpath(genpath('../')); % add entire project path, incl. general code used across analyses

%% Parameters 
% params for analysis type
signal_type = 'envl.dat'; % 'muae'
runV1 = false;
runArch = false;

% params for sparseNoiseExpByChan
zcutoffs = inf; %[8 10 20];
RFtwin = [0 0.15];
saveRaw = true;
plotFlag = true;
saveFigs = true;
saveMats = true;

% params for getAllunits
corr_p      = 0.001;         % 0.001; corr_p represents the wilcoxon ransum correlation for between/withing segments
qi          = 0.05;          % 0.05; Philipps quality index
cluster_qi  = 3;             % 3; cluster quality index of likely single- vs multiunit
r           = 0.65;          % 0.7; regression index
type        = 3;

% Save directory
if runV1
    dirOfDirs = '/Volumes/busse_lab/users/sinem/SparseNoise_V1';    
else    
%     dirOfDirs = '/Volumes/lab/users/miro/data/RF';
    dirOfDirs = '/Volumes/lab/users/yannik/projects/ret2dlgn/analyses/dlgn_rf_depth/';    
end

if runArch
    dirOfDirs = [dirOfDirs 'arch'];
end

%% Select Mice
% GET KEYS USED FOR UNITS IN OTHER ANALYSES OF ROMAN_ROSON_ET_AL_2018
%myMice.mouse_id = 'BL6_0190';

myMice = fetch(data.Mice(getAllUnits(corr_p, qi, cluster_qi, r, type)),'mouse_id');
if ~isdir(dirOfDirs)
    mkdir(dirOfDirs);
end

%myMice = {'Ntsr1-Cre_0082', 'Ntsr1-Cre_0080'}; % chr2patch
%myMice = {'Ntsr1-Cre_0088', 'Ntsr1-Cre_0085'}; % archpatch

%% Analyze files
for zcutoff = zcutoffs
    if zcutoff > 8 && saveRaw
        saveRaw = false;
    end
    
    for midx = 1:numel(myMice)
        myMouseID = myMice(midx).mouse_id;
        myAnimal = fetch(data.Mice(sprintf('mouse_id like "%s"', myMouseID)));
        
        if runV1
            if runArch
                v1Series = fetch(data.Series(myAnimal,'series_depth < 1800'));
            else
                v1Series = fetch(data.Series(myAnimal,'series_depth < 1800') & data.ClusterInfo);
            end
            lgnSeries = v1Series;
        else
            if runArch
                lgnSeries = fetch(data.Series(myAnimal,'series_depth > 1800'));
            else
                lgnSeries = fetch(data.Series(myAnimal,'series_depth > 1800') & data.ClusterInfo & data.Series('series_date > "2014-03-14"'));
            end
        end

        %%% %%%% %%%% %%%% %%%% 
        %%% FOR DEBUGGING %%%%        
        % lgnSeries = lgnSeries(ismember([lgnSeries.series_num] , 6));        
        %%% FOR DEBUGGING  OR SELECTING A SERIES %%%%
        %%% %%%% %%%% %%%% %%%% 
        
        if runArch
            myExps = fetch(data.Experiments('exp_name like "%sparse%"', lgnSeries));
        else
            myExps = fetch(data.Experiments('exp_name like "%sparse%"', lgnSeries) & data.Units);
        end
        
        if saveFigs
            figDir = fullfile(dirOfDirs, 'results', sprintf('sparseNoiseFigs_%s', upper(signal_type)));
            if ~isdir(figDir)
                mkdir(figDir)
            end
        else
            figDir = [];
        end
        
        if saveMats
            matDir = fullfile(dirOfDirs, 'data', sprintf('sparseNoiseMats_%s', upper(signal_type)));
            if ~isdir(matDir)
                mkdir(matDir)
            end
            
        else
            matDir = [];
        end

        % animal80,ser11exp11 - should be 9 or 10?;
        
        % find out which experiments to combine
        skeys = lgnSeries;
        nser = numel(skeys);
        uniqueExps = cell(1, nser);

        % needed to rerun Ntrs1-Cre_0080s07 bc DB thought that there was
        % opto in some spare noise experiments. 
        % reRunIdx = find([skeys.series_num] == 7);
        % reRunIdx = find([skeys.series_num] == 9);
        
        for iseries = 1:nser
            % all exps in series
            ekeys = myExps([myExps.mouse_counter] == skeys(iseries).mouse_counter &...
                [myExps.series_num] == skeys(iseries).series_num);
            % get the unique experiments (combine exps that are the same)
            
            try
                uniqueExps{iseries} = findSeriesExp(myMouseID, sprintf('%s', '%sparse%'), skeys(iseries).series_num);
            catch
                warning('findSeriesExp:mouse_id: %s and series: %d',myMouseID,skeys(iseries).series_num)
                continue
            end
            % number of unique exps
            nexps = numel(uniqueExps{iseries});
            
            % for each  unique exp
            for ie = 1:nexps
                %  keys of exps
                ikey = ekeys(ismember([ekeys.exp_num], [uniqueExps{iseries}(ie).exp_num]));
                
                % DEBUGGING: test keys
%                 key(1).mouse_counter = 132;
%                 key(1).series_num = 12;
%                 key(1).exp_num = 1;
%                 key(2).mouse_counter = 132;
%                 key(2).series_num = 12;
%                 key(2).exp_num = 6;
                
%                 key(1).mouse_counter = 58;
%                 key(1).series_num = 6;
%                 key(1).exp_num = 1;
                
%                 key(2).mouse_counter = 58;
%                 key(2).series_num = 6;
%                 key(2).exp_num = 3;
                
                % run sparsenoise code               
                try
%                     sparseNoiseExpByChan(ikey, RFtwin, zcutoff, signal_type, plotFlag, figDir, matDir,saveRaw)                
                    sparseNoiseExpByChan(ikey, RFtwin, zcutoff, signal_type, plotFlag, figDir, matDir, saveRaw)
                catch
                    warning('mouse_id %s; series %d; exp_num: %d',myMouseID,ikey.series_num,ikey.exp_num)
                end
                
                % pause
                close all;
                
                % also do it for each experiment separately but only if it was
                % combined (best here for saving but better outside when
                % not saving for comparison)
                if numel(ikey) > 1
                    for iie = 1:numel(ikey)
                        eMatDir = fullfile(matDir, 'eachExperiment');
                        eFigDir = fullfile(figDir, 'eachExperiment');
                        if ~isdir(eMatDir); mkdir(eMatDir);end
                        if ~isdir(eFigDir); mkdir(eFigDir);end
                        try
                            sparseNoiseExpByChan(ikey(iie), RFtwin, zcutoff, signal_type, plotFlag, eFigDir, eMatDir,[]);
                        catch
                            % continue
                        end
                        %pause
                        close all;
                    end
                end
            end
        end

        %pause
        %sound(achtung)
    end
end

return
