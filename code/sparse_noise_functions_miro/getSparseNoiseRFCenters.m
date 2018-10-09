
function getSparseNoiseRFCenters(mouse_id, mySeries, pickNewRfs, select_z_cutoff, alwaysPlotIntRF, muaDir, figDir,rfDir, tightPlotFlag, figPositions)
% load sparseNoise Mat created in makeSparseNoiseMats.m and lets you either
% pick RF centers manually OR if you already did/RF centers exist in some
% rfDir, then it lets you plot them;
% NEEDS setFigShit to be in the path
%
% getSparseNoiseRFCenters(mouse_id, lgnSeries, pickNewRfs, select_z_cutoff, alwaysPlotIntRF, muaDir, figDir,rfDir, tightPlotFlag)
%
% ONLY DOES ONE MOUSE AT A TIME
%
% DEFAULTS:
% mySeries - lgnSeries: all series with series depth > 1800  & data.ClusterInfo
% pickNewRfs - true
% select_z_cutoff - true
% alwaysPlotIntRf - true
% muaDir - '/Volumes/busse_lab/users/sinem/SparseNoise_LGN/sparseNoiseMats_MUAE';
% figDir - '/Volumes/busse_lab/users/sinem/SparseNoise_LGN/sparseNoiseFigs_MUAE';
% rfDir - fullfile(muaDir, sprintf('manualRFs_%s', datestr(datetime('now'),'yyyymmdd')))
%       - e.g.   fullfile(muaDir, 'manualRFs_20160127')
%       - note that if you try to run it multiple times on the same day and
%       choose not to overwrite when prompted, it will make a new directory
%       fullfile(muaDir, 'manualRFs_20160127(1)')
% tightPlotFlag - false (whether to use subplottight instead of subplot)
% figPositions  - cell array that has {rfPos and zPos}, useful when comparing picked RFs from different directories
%               - by default: rfPos = [0 30 60] and zPos = [33, 16, 0];
%               - rfPos is the x-pos (on screen) of combined, on and off figures
%               - zPos is the y-pos (on screen of the different z scores when select_z_cutoff=T and pickNewRfs=T)
%               - if pickNewRfs=F, zPos(2) is where the RFs are plotted
%               - if pickNewRfs=F & alwaysPlotIntRF=T, zPos(1) is where the IntRFs are plotted
%
%
% see also: chooseAndVisualizeSparseNoiseRFs, makeSparseNoiseMats, sparseNoiseExpByChan


setFigShit;

% MOUSE/SERIES TO WORK ON
myAnimal = fetch(data.Mice(sprintf('mouse_id like "%s"', mouse_id )));
if nargin < 2 || isempty(mySeries)
    mySeries = fetch(data.Series(myAnimal,'series_depth > 1800') & data.ClusterInfo);
end
% COMPUTE NEW RFs? otherwise will just plot them
if nargin < 3 || isempty(pickNewRfs)
    pickNewRfs = true;
end
% SELECT Z-CUT OFF WHILE PICKING RFS? If not, will always use 8 besides 2
% special cases (see code)
if nargin < 4 || isempty(select_z_cutoff)
    select_z_cutoff = true;
end
% if just plotting the RFs, also plot the IntRFByCh if another rfEstimate
% was used?
if nargin < 5 || isempty(alwaysPlotIntRF)
    alwaysPlotIntRF = true;
end

% muaDIR - where the sparseNoiseMat files reside and where the rfDir will
% be placed
if nargin < 6 || isempty(muaDir)
    muaDir = '/Volumes/busse_lab/users/sinem/SparseNoise_LGN/sparseNoiseMats_MUAE';
end
% figDIR - where the figures were saved during running makeSparseNoiseMats
% (and call to sparseNoiseExpByCh)
if nargin < 7 || isempty(figDir)
    figDir = '/Volumes/busse_lab/users/sinem/SparseNoise_LGN/sparseNoiseFigs_MUAE';
end

% IF pickNewRfs IS FALSE SPECIFY rfDIR
if (nargin < 8 || isempty(rfDir))
    if ~pickNewRfs
        error('getSparseNoiseRFCenters:noDirNoPlot', 'plotting of RF centers requires rfDir to be specified');
    else
        rfDir = fullfile(muaDir, sprintf('manualRFs_%s', datestr(datetime('now'),'yyyymmdd'))); %#ok<*UNRCH>
    end
end

if nargin < 9 || isempty(tightPlotFlag)
    tightPlotFlag = false;
end

% DEFAULT FIGURE POSITIONS
if nargin < 10 || isempty(figPositions)
    rfPos = [0 30 60];
    zPos = [33, 16, 0];
else
    if isempty(figPositions{1})
        rfPos = [0 30 60];
    else
        rfPos = figPositions{1};
    end
    
    if isempty(figPositions{2})
        zPos = [33, 16, 0];
    else
        zPos = figPositions{2};
    end
    
end

%% DEFAULT RF TYPES, AND DEFAULT RF ESTIMATE
% z cut offs to cycle through (only relevant when pickNewRfs=T)
zCos = [8, 10, 20];
% RF TYPES CAN BE '' (Combined), 'ON' OR 'OFF'
rfType = {'', 'on','off'};
% rfEstimate should always be IntRFByCh, however, there is now an option to use
% the diff if all else fails but only if select z cutoff is set to true
% - used for Ntsr1-Cre_0082 series 4
rfEstimate = 'IntRFByCh'; % SHOULD NOT REALLY CHANGE!



%% RF DIR ASSIGNMENT
if pickNewRfs
    if exist(rfDir,'dir')
        warning('Directory %s already exists', rfDir);
        dirFiles = dir(fullfile(rfDir,sprintf('*manualRFcenters.mat')));
        % check to see if files already exist in the directory for each series
        serFiles = arrayfun(@(s) sprintf('%s_%d_manualRFcenters.mat', mouse_id, s.series_num), mySeries, 'uniformoutput', 0);
        
        seriesFileAlreadyPresent = arrayfun(@(serfile) any(strcmpi(serfile, {dirFiles.name})), serFiles);
        
        % ask to overwrite if they do
        if any(seriesFileAlreadyPresent)
            warning('Files for %d out of %d specified series already exist', nnz(seriesFileAlreadyPresent), numel(mySeries));
            userIn = input('Shall i create a new directory? (Ctrl-C to exit) \n y/n? ', 's');
            newdir = strcmpi(userIn, 'y');
            if newdir
                [~,dir_version] = regexp(rfDir,'\((\d+)\)$', 'match', 'once', 'tokens');
                if isempty(dir_version)
                    dir_version = 0;
                    rfDir = sprintf('%s(%d)', rfDir, dir_version+1);
                else
                    dir_version = str2double(dir_version);
                    rfDir = regexprep(rfDir, '\d+(?=\)$)', num2str(dir_version+1));
                end
            elseif ~all(seriesFileAlreadyPresent)
                overwrite_prompt = sprintf('Do you want to continue processing the next series? (if no, I overwrite all %d files) y/n', nnz(seriesFileAlreadyPresent));
                userIn = input(overwrite_prompt,'s');
                if strcmpi(userIn, 'y')
                    mySeries= mySeries(~seriesFileAlreadyPresent);
                end
            end % new dir
        end % rfDir contains files for specified series
    end % rfDir already exists
    
    if ~exist(rfDir,'dir')
        mkdir(rfDir);
    end
    
elseif ~exist(rfDir, 'dir')
    error('Specified directory %s does not exist', rfDir);
end % pick new RFs

%% COMPUTE NEW RF CENTERS
if pickNewRfs
    for iser = 1:numel(mySeries) %1:numel(lgnSeries)
        series_num = mySeries(iser).series_num;
        % zcutoff
        if ~ select_z_cutoff
            if      strcmp(mouse_id, 'Ntsr1-Cre_0082') && series_num == 7 || ...
                    strcmp(mouse_id, 'Ntsr1-Cre_0080') && series_num == 3
                z_cutoff = 10;
                
            else
                z_cutoff = 8;
            end
            tmp = dir(fullfile(figDir, sprintf('z_co%d', z_cutoff), sprintf('*Map_*%s_%02d*', mouse_id,series_num)));
        else
            figDir = sprintf('/Volumes/busse_lab/users/sinem/SparseNoise_LGN/sparseNoiseFigs_MUAE');
            tmp = dir(fullfile(figDir, 'z_co8', sprintf('*Map_*%s_%02d*', mouse_id,series_num)));
        end
        
        [~, tmp1] = regexp({tmp.name}, '_([\S]+).fig', 'match', 'tokens', 'once');
        basenames = [tmp1{:}];
        nBase = numel(basenames);
        clear tmp1 tmp
        
        
        serRFs = struct('mouse_id', repmat({mouse_id}, 1, nBase), ...
            'series_num', num2cell(repmat(series_num, 1, nBase)), ...
            'exps',cell(1, nBase ), 'z_cutoff',cell(1, nBase ), 'rfType', cell(1,nBase),...
            'RFcenters', cell(1, nBase ), 'path2map', cell(1,nBase), 'rfEstimate', repmat({rfEstimate}, 1, nBase));
        
        for ibase = 1:nBase
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% DEFAULT RF ESTIMATE TYPE %%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~strcmpi(rfEstimate, 'IntRFByCh');
                rfEstimate = 'IntRFByCh'; % SHOULD NOT REALLY CHANGE!
            end
            
            % base name (mouseid, series and exps)
            baseStr = basenames{ibase};
            % assign experiments to serRFs struct
            [~, tmp] = regexp(baseStr, '\(([\S]+)\)', 'match', 'tokens', 'once');
            serRFs(ibase).exps = str2num(tmp{:}); %#ok<ST2NM>
            
            if select_z_cutoff
                for iz = 1:numel(zCos)
                    for itype = 1:numel(rfType)
                        % load figure from figDir
                        ifname = fullfile(figDir, sprintf('z_co%d', zCos(iz)), sprintf('sparseNoiseMap%s_%s.fig', rfType{itype}, baseStr));
                        hf = open(ifname);
                        % new name of figure that includes zco
                        iFigStr = sprintf('z_co%d \t %s', zCos(iz), hf.Name);
                        % position and rename figure
                        set(hf, figPars, 'position', [rfPos(itype) zPos(iz), 20,15], 'Name', iFigStr);
                        nm = ceil(sqrt(numel(hf.Children)));
                        make_subplottight(hf, nm,nm)
                        set(get(hf,'Children'), 'xtick', [], 'ytick', [])
                        clear hf
                    end
                end
                z_cutoff = [];
                while isempty(z_cutoff)
                    zprompt = sprintf('Please specify the best z_cutoff. Pick one of scalar inputs %s. \n (or `diff`,`zbych`,`diff_zbych` to change rf estimate - or `reset` to choose from original) \n', strjoin(cellfun(@num2str, (num2cell(zCos)), 'uniformoutput',0), ','));
                    zin= input(zprompt);
                    z_cutoff = zCos(ismember(zCos, zin));
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%% CASE WHERE DIFF RECEPTIVE FIELD ESTIMATE IS USED IS USED! %%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if strcmpi(zin, 'diff') || strcmpi(zin, 'zbych') || strcmpi(zin, 'diff_zbych')
                        warning('You are now using a different method of estimating the RF... this may take a while...')
                        if strcmpi(zin, 'diff')
                            rfEstimate = 'IntdRF';
                        elseif strcmpi(zin, 'zbych')
                            rfEstimate = 'IntRFByCh_zbych';
                        elseif strcmpi(zin, 'diff_zbych')
                            rfEstimate = 'IntdRF_zbych';
                        end
                        
                        for iz = 1:numel(zCos)
                            diffSearch = fullfile(muaDir, sprintf('z_co%d', zCos(iz)),[baseStr '*Info*']);
                            diffFile = dir(diffSearch);
                            if numel(diffFile)~=1
                                error('why cant i load the the right rfVars for the diff case now? multiple files!')
                            end
                            diffVars = load(fullfile(muaDir, sprintf('z_co%d', zCos(iz)), diffFile.name));
                            for itype = 1:numel(rfType)
                                diffRF = eval(sprintf('diffVars.%s%s', rfType{itype}, rfEstimate));
                                figure(figPars, 'Name', sprintf('%s_%sRF - %s RF!', baseStr , rfType{itype}, upper(zin)), 'position', [rfPos(itype)+15 zPos(iz), 20,15]);
                                plotRFbyChStruct(diffRF, diffVars.sRegionX, diffVars.sRegionY,diffVars.intRegionX,diffVars.intRegionY,tightPlotFlag);
                            end % itype
                        end % iz
                    elseif strcmpi(zin, 'reset')
                        rfEstimate = 'IntRFByCh';
                    end  % case where diff rfEstimate is used
                    %%% END CASE WEHRE DIFF RF ESTIMATE IS USED %%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                end % while isempty( z-cutoff)
                
                % assign rfEstimate to serRFs struct
                serRFs(ibase).rfEstimate = rfEstimate;
            end % if select zcutoff
            
            
            % assign z_cutoff to serRFs structure
            serRFs(ibase).z_cutoff = z_cutoff;
            % get the directory where the matfile resides
            matDir = fullfile(muaDir, sprintf('z_co%d', z_cutoff));
            % list the files that match the particular basename
            flist = dir(fullfile(matDir, [baseStr '*']));
            % there should only be 1
            if numel(flist) > 1
                flist = flist(cellfun(@isempty, strfind({flist.name},'RFcenters')));
                if numel(flist) > 1
                    error('multiple mat files for ls %s?',fullfile(matDir, [baseStr '*']));
                end
            end
            fname = flist.name;
            
            close all;
            % load file
            rfVars = load(fullfile(matDir,fname));
            % assign path to file to serRFs struct
            serRFs(ibase).path2map = fullfile(matDir,fname);
            % number of channels
            nch = numel(rfVars.IntRFByCh);
            % initialize
            RFcenters = nan(nch,2);
            for itype = 1:numel(rfType)
                h(itype) = figure(figPars, 'Name', sprintf('%s_%sRF - SELECT %d', ...
                    baseStr , rfType{itype}, itype), 'position', [rfPos(itype) zPos(2), 30,25]); %#ok<AGROW>
                
                thisRF = eval(sprintf('rfVars.%s%s', rfType{itype}, rfEstimate));
                
                % plot and get axes
                AH(itype).a = plotRFbyChStruct(thisRF, rfVars.sRegionX, rfVars.sRegionY,rfVars.intRegionX,rfVars.intRegionY,tightPlotFlag); %#ok<AGROW>
                
            end
            % ask which RF type to use
            figLog = false(1,numel(rfType));
            myType = input('Which RF type is the best? Press 1 for combined, 2 for on, 3 for off');
            % close all other figures and assign the best RF figure as the current figure
            figLog(myType) = true;
            set(h(figLog), 'position', [rfPos(3) zPos(2), 30,25])
            ah = AH(figLog).a;
            close(h(~figLog));
            h = h(figLog);
            
            % assign best RF type to serRFs struct
            serRFs(ibase).rfType = rfType{myType};
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% START USER TASK OF FIND THE CENTERS %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            userIn = input('Now its time to find RF centers. For each channel that has an RF, click on the center. Press "c" to begin. If no RFs are present, press "d"', 's');
            if strcmpi(userIn,'d')
                userIn = input('Press "d" again to confirm or "c" to continue \n', 's');
            end
            
            while strcmpi(userIn,'c')
                figure(h);
                ich = [];
                while isempty(ich)
                    [~] = input('click on channel and press "c" when ready.\n','s');
                    ich = find(ah==gca);
                    chanLabel  = get(get(ah(ich),'Title'), 'String');
                    xhairs = findobj(ah(ich), 'color','r');
                    if ~isempty(xhairs)
                        delete(xhairs);
                    end
                end
                fprintf('Select rf center on CH%s and press enter. Only the last coordinates will be used. \n', chanLabel);
                h2 = figure(figPars, 'Name', sprintf('CH%s',chanLabel),'position', [rfPos(2) zPos(2), 30,25]);
                a2 = axes; copyobj(get(ah(ich), 'Children'), a2);
                colormap gray; axis square; set(a2, 'XTick', [], 'YTick', [], 'YDir', 'Normal');
                [x,y] = getpts(a2);
                % so that you can not take a closer look and decide no RF
                if isempty(x) || isempty(y)
                    x=nan; y=nan;
                end
                RFcenters(ich, :) = [x(end) y(end)];
                close(h2);
                axes(ah(ich)); hold on; %#ok<LAXES>
                plot(ah(ich), x(end), y(end), 'r+', 'markersize', 5);
                userIn = input('Press "c" to continue to the next channel "d" if you are done \n', 's');
                
                if strcmpi(userIn,'d')
                    userIn = input('Press "d" again to confirm or "c" to continue \n', 's');
                end
                
            end
            
            % assign RF centers to serRFs struct
            serRFs(ibase).RFcenters = RFcenters;
            clear RFcenters;
        end % ibase
        fprintf('saving serRFs struct for %s series %d \n', mouse_id, series_num);
        save(fullfile(rfDir, sprintf('%s_%d_manualRFcenters', mouse_id, series_num)), 'serRFs');
    end % iser
    
    %% PLOT PREVIOUSLY FOUND RF CENTERS
else
    for iser = 1:numel(mySeries) %1:numel(lgnSeries)
        series_num = mySeries(iser).series_num;
        % load serRFs
        tmp = load(fullfile(rfDir, sprintf('%s_%d_manualRFcenters', mouse_id, series_num)));
        serRFs = tmp.serRFs; clear tmp;
        for iexp = 1:numel(serRFs)
            % load the RF information for the experiment
            rfVars = (load(serRFs(iexp).path2map));
            baseStr = regexp(serRFs(iexp).path2map, '[^/]+(?=_RFInfo)', 'match','once');
            
            for itype = 1:numel(rfType)
                % MARKER IS RED FOR THE RF THAT WAS USED AND GREEN OTHERWISE
                if strcmp(serRFs(iexp).rfType,rfType{itype})
                    imarker = 'r+';
                else
                    imarker = 'go';
                end
                
                % PLOT RF THAT WAS USED
                figure(figPars, 'Name', sprintf('z_co=%d \t %s_%sRF - USED %s', ...
                    serRFs(iexp).z_cutoff, baseStr , rfType{itype}, serRFs(iexp).rfEstimate), ...
                    'position', [rfPos(itype) zPos(2), 20,15]);
                thisRF = eval(sprintf('rfVars.%s%s', rfType{itype}, serRFs(iexp).rfEstimate));
                % plot RFs and get axes
                ah = plotRFbyChStruct(thisRF, rfVars.sRegionX, rfVars.sRegionY,rfVars.intRegionX,rfVars.intRegionY,tightPlotFlag);
                set(ah, 'NextPlot', 'add')
                arrayfun(@(ich) plot(ah(ich), serRFs(iexp).RFcenters(ich,1),serRFs(iexp).RFcenters(ich,2), imarker, 'markersize', 5), 1:size(serRFs(iexp).RFcenters,1));
                
                clear ah
                % IF USED RF is not IntRFByCh, then plot IntRFByCh as well
                if ~strcmpi(serRFs(iexp).rfEstimate, 'IntRFByCh') && alwaysPlotIntRF
                    imarker = 'go';
                    figure(figPars, 'Name', sprintf('z_co=%d \t %s_%sRF - IntRFByCh', ...
                        serRFs(iexp).z_cutoff, baseStr , rfType{itype}), ...
                        'position', [rfPos(itype) zPos(1), 20,15]);
                    thisRF = eval(sprintf('rfVars.%s%s', rfType{itype}, 'IntRFByCh'));
                    % plot RFs and get axes
                    ah = plotRFbyChStruct(thisRF, rfVars.sRegionX, rfVars.sRegionY,rfVars.intRegionX,rfVars.intRegionY,tightPlotFlag);
                    set(ah, 'NextPlot', 'add')
                    arrayfun(@(ich) plot(ah(ich), serRFs(iexp).RFcenters(ich,1),serRFs(iexp).RFcenters(ich,2), imarker, 'markersize', 5), 1:size(serRFs(iexp).RFcenters,1));
                end % if IntRFByCh was not used
            end % for rfType
            drawnow
        end % for iexp
    end % for iser
end % if pickNewRfs


