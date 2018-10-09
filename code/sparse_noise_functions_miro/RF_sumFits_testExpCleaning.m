% takes chol decomp fits and ...
% Xf = on fits
% Yf = off fits
% Zf = on + off fits
% P(x+y|Xf,Yf, Zf)= ? Don't do this.
% weighted sum - mixture of gaussains
% aXr + bYr + cZr, where a,b, and c are the weights (either equal or depends on R2)
% dont sum distribtions, instead samples from them and then kmeans k=1
% how do you do the weighting? take a total of N samples, and weight contribution


% rsq co (to determine good fit)
close all; clear all; startup;
addpath ~/Dropbox/clean/locomotion_layer/animals/agneReAnalysis/sparseNoiseFunctions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% options  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
runV1 = false;
onlyCleanedExps = false;
plotSeriesStruct = true;
makeSeriesStruct = false;
makeSumStruct = false;
myRsqCutoffs = [0.45 0.5 0.60] ;% for V1 only run 0.60  %[0.40 0.45 0.50 0.60];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% use chol decomp fitting method? (fitGuass2D_titled is much better)
% thisFit = 'fitGauss2D_chol';
thisFit = 'fitGauss2D_tilted';

% directory of original fits (should really be a table...)
if runV1
    %muaDir = '/Users/yoknapatawpha/Desktop/SparseNoise_V1/sparseNoiseMats_MUAE';
    muaDir = '/Volumes/busse_lab/users/sinem/SparseNoise_V1/sparseNoiseMats_MUAE';
else
    %muaDir = '/Users/yoknapatawpha/Desktop/SparseNoise_LGN/sparseNoiseMats_MUAE';
    muaDir = '/Volumes/busse_lab/users/sinem/SparseNoise_LGN/sparseNoiseMats_MUAE';
end
% directory where invidivual fits per experiment type reside
fitDir = fullfile(muaDir, thisFit);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% ONLY PLOT SERIES STRUCT ? %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotSeriesStruct = true; % does this with different Rsq options, see below
% for visualization
setFigShit
zPos = [33, 16, 0];
rfPos = [0 30 60];
subplotsAreTight = true;
scrnfigw = 24;
scrnfigh = 16;
setFigShit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% MAKE SUM STRUCT ? %%%%%%%%%%%%%%%
% note: seriesStructs should be computed first
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%makeSumStruct = false; % now definied above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% MAKE SERIES STRUCT ? %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% struct where experiments are concatenated (not summed yet)
% find best zco and resave structures for this type of fitting
% makeSeriesStruct = false;% now definied above

% all RFtypes
rfType = {'', 'on', 'off'};
rfTypeFullString = cellfun(@(type) strcat(type, 'IntRFByCh'), rfType, 'uniformoutput', 0);
nrf = numel(rfType);

% all z-cutoffs (best one will be picked for each experiment type)
z_cos = [8,10,20];
nzcos = numel(z_cos);

% get all basenames
tmpDir = fullfile(muaDir, sprintf('z_co%d',8));
tmp = dir(fullfile(tmpDir, '*_RFInfo*'));
[~, tmp1] = regexp({tmp.name}, '([\S]+)\_RFInfo.mat', 'match', 'tokens', 'once');
basenames = [tmp1{:}];

% get all animals and series aligned to basenames
allMice = regexp(basenames, 'Ntsr1-Cre_[0-9]{4}', 'match', 'once');
[tmp,~] = regexp(basenames, '_([0-9]){2}_', 'tokens', 'match', 'once');
allSeriesString = [tmp{:}];
allSeries = cellfun(@str2num, allSeriesString);
unqMice = unique(allMice);



if makeSumStruct || makeSeriesStruct
    for rsq_co = myRsqCutoffs;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%  FOR SUM STRUCT   %%%%%%%%%%%%%%%
        % note: seriesStructs should be computed first
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % directory for summed gaussian structs (depends on rsq)
        rfCdir = fullfile(muaDir, sprintf('summedFits_rsq%d', rsq_co*100));
        if ~isdir(rfCdir)
            mkdir(rfCdir);
        end
        
        sumDir = fullfile(rfCdir, thisFit);
        if ~isdir(sumDir)
            mkdir(sumDir);
        end
        
        if onlyCleanedExps
            sumSerDir = fullfile(sumDir,'cleanExpFits');
        else
            sumSerDir = fullfile(sumDir,'seriesFits'); %#ok<*UNRCH>
        end
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%   FOR SERIES STRUCT  %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % new structs will be saved under:
        % directory for series structs
        serDir = fullfile(sumDir, 'seriesInfo');
        % TODO : for different rsqs, new series structs will be saved which is a
        % bit of a waste of space ...
        
        
        
        
        %% makeSeriesStruct
        if makeSeriesStruct || ~isdir(serDir)
            rng(1); % for reproducibility
            
            fprintf('\t Making seriesStruct \n');
            
            if ~isdir(serDir)
                mkdir(serDir);
            end
            
            % all zcos directories for the old structs
            zFitDirs = arrayfun(@(z) fullfile(fitDir, sprintf('z_co%d',z)),...
                z_cos, 'uniformoutput', 0);
            
            for imouse = 1:numel(unqMice)
                
                mouse_id = unqMice{imouse};
                % aligned to basenames
                mouseLog = ~cellfun(@isempty, strfind(basenames, mouse_id));
                mouseSeries = unique(allSeries(mouseLog));
                
                % mouse counter
                mc = fetch1(data.Mice(sprintf('mouse_id like "%s"', mouse_id)), 'mouse_counter');
                
                for iser = 1:numel(mouseSeries)
                    
                    thisSer = mouseSeries(iser);
                    
                    % series key
                    skey = struct('mouse_counter', mc, 'series_num', thisSer);
                    % indices into basenames
                    serIndices = find(mouseLog & allSeries==thisSer);
                    % number of exp types in series
                    nbase = numel(serIndices);
                    
                    
                    % experiment types
                    [~,baseExps] = arrayfun(@(i) ....
                        regexp(basenames{i}, '\(([0-9,\,]+)\)', 'match', 'tokens', 'once'),...
                        serIndices, 'uniformoutput', 0);
                    baseExps = cellfun(@(ibase) cellfun(@str2num,strsplit(ibase, ',')), ...
                        [baseExps{:}], 'uniformoutput',0);
                    
                    
                    % seriesStruct is organized zco x nbaseExps
                    tmpS = struct('mouse_id', mouse_id, 'series_num', thisSer,'exps', {[]}, ...
                        'z_cutoff', {[]},'fitType',thisFit,'expInfo', {[]}, ...
                        rfTypeFullString{1},{[]}, rfTypeFullString{2},{[]}, rfTypeFullString{3},{[]});
                    seriesStruct = repmat(tmpS, nzcos, nbase);
                    % add in exps and z_cutoff
                    tmpe = repmat(baseExps,nzcos,1); tmpz = num2cell(repmat(z_cos',1,nbase));
                    [seriesStruct.exps] = deal(tmpe{:}); clear tmpe;
                    [seriesStruct.z_cutoff] = deal(tmpz{:}); clear tmpz;
                    
                    
                    % sparsenoise experiment keys
                    ekeys = fetch(data.Experiments & data.SparseNoiseTuning(skey));
                    
                    
                    % get monitor positions (use first experiment of each group)
                    monitorExps = cellfun(@(s) s(1), baseExps);
                    monitorExps = fetch(data.Experiments(skey, sprintf('exp_num in (%s)',strjoin(arrayfun(@(e) num2str(e), monitorExps, 'uniformoutput', 0),','))));
                    [mon_xx,mon_yy,stim_xc,stim_yc, monitor_xy] = arrayfun(@(exp) ....
                        calculateStimulusSpace(exp), monitorExps, 'uniformoutput', 0);
                    % get the xcoords and ycoords of the rfmaps
                    % load one z_co to get x,y coords
                    zDir1 = fullfile(muaDir, sprintf('z_co%d',z_cos(1)));
                    tmp = arrayfun(@(idx) load(fullfile(zDir1, sprintf('%s_RFInfo.mat',basenames{idx}))),serIndices);
                    
                    % put it all monitor info in a structure
                    monitorStruct = struct('xcoords',{tmp.intRegionX}','ycoords',{tmp.intRegionY}',...
                        'stim_xc', stim_xc, 'stim_yc', stim_yc, ...
                        'mon_ang', cellfun(@(xy) xy(1), monitor_xy, 'uniformoutput', 0), ...
                        'mon_elv', cellfun(@(xy) xy(2), monitor_xy, 'uniformoutput', 0), ...
                        'mon_xx', mon_xx, 'mon_yy', mon_yy); clear tmp;
                    monitorStruct = repmat(monitorStruct(:), 1,nzcos)';
                    
                    % assign them to seriesStruct
                    tmp = arrayfun(@(i) i, monitorStruct, 'uniformoutput', 0);
                    [seriesStruct.expInfo] = deal(tmp{:}); clear tmp;
                    
                    % put in the fit data and get summary info (e.g. best z_co)
                    for ibase = 1:nbase
                        
                        % for this experiment type, load fits for all z_cos
                        zfiles = cellfun(@(zdir) fullfile(zdir, [basenames{serIndices(ibase)},'_fits.mat']), ...
                            zFitDirs, 'uniformoutput', 0);
                        
                        tmp = cellfun(@(file) load(file), zfiles);
                        zFits = cat(1,tmp.myFit);
                        
                        % number of channels
                        nch = size(zFits,2);
                        
                        % assign zFits to seriesStruct
                        for irf = 1:nrf
                            rfString = rfTypeFullString{irf};
                            for iz = 1:nzcos
                                seriesStruct(iz,ibase).(rfString) = [zFits(iz,:).(rfString)];
                            end
                        end
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%% GET SUMMARY INFO! (e.g. about best z_co)%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % check if nan across all RF types
                        nanZlog = any(cell2mat(cellfun(@(rf) any(arrayfun(@(z) ...
                            isnan(z.(rf).x0), zFits),2), rfTypeFullString, 'uniformoutput', 0)),2);
                        
                        % for indexing the base
                        base_zcos = z_cos(~nanZlog);
                        single_zco = nnz(~nanZlog) == 1;
                        
                        % make a rsqMat [nZ, nRF, nCh]
                        %rsqMat = arrayfun(@(z) cellfun(@(RF) zFits(1).(RF).r2, rfTypeFull), zFits, 'uniformoutput',0);
                        tmp = arrayfun(@(ch) cell2mat(arrayfun(@(z) cellfun(@(RF) z.(RF).r2, rfTypeFullString), zFits(:,ch), 'uniformoutput',0)), 1:nch , 'uniformoutput',0);
                        rsqMat = cat(3,tmp{:});clear tmp;
                        % for each z_co, figure out how many channels > rsq_co
                        goodRFcounts_all = squeeze(sum(rsqMat > rsq_co,2));
                        goodRFcounts = squeeze(sum(rsqMat(~nanZlog,:,:)  > rsq_co,2));
                        
                        % number of channels that pass
                        nChPass_all = sum(goodRFcounts_all > 0, 2);
                        nChPass = sum(goodRFcounts > 0, (~single_zco) + 1);
                        
                        % if no channels passed then set everything to 0
                        if ~any(nChPass)
                            [z_bestRsq,z_mostRF,z_mostCh] = deal(0);
                            % if only one base_zcos is left, then no zco selection needed
                        elseif single_zco
                            z_bestRsq = base_zcos;
                            z_mostCh = base_zcos;
                            z_mostRF = base_zcos;
                            %tmp = arrayfun(@(ch) cell2mat(arrayfun(@(z) cellfun(@(RF) z.(RF).r2, rfTypeFullString), zFits(:,ch), 'uniformoutput',0)), 1:nch , 'uniformoutput',0);
                        else
                            z_mostCh = base_zcos(nChPass == max(nChPass));
                            % z_co that has the most chans > rsq_co across RFtypes
                            nRFPass = sum(goodRFcounts,2);
                            z_mostRF = base_zcos(nRFPass  == max(nRFPass));
                            
                            % z-score that has the best fit (regardless of rf type) for each ch
                            maxR2 = squeeze(max(rsqMat(~nanZlog,:,:),[],2)); % best fit (across RF type) for each zscore
                            [~,tmpi] = max(maxR2); % tmpi is an index into base_zcos
                            z_bestRsq = mode(nonzeros(base_zcos(tmpi) .* any(goodRFcounts)));clear tmpi;
                        end % isscaler(base_zcos)
                        
                        
                        % structure for fields to be added at the level of seriesStructs
                        tmpS = struct('anyNan', num2cell(nanZlog'), 'nChPass', num2cell(nChPass_all'), ...
                            'mostCh', num2cell(ismember(z_cos, z_mostCh)), ... % logical
                            'mostRF', num2cell(ismember(z_cos, z_mostRF)), ... % logical
                            'bestRsq', num2cell(ismember(z_cos, z_bestRsq)) ... % logical
                            )' ;
                        
                        
                        % does the one with nan have more channels?
                        % know it happens for  'Ntsr1-Cre_0080_05_(1)'
                        nanHasMoreCh = nChPass_all(nanZlog) > nChPass(base_zcos == z_mostCh(1));
                        if isempty(nanHasMoreCh); nanHasMoreCh = false; end;
                        
                        % structure for fields to be added at the level of expInfo (include for all Zs)
                        tmpI = struct('z_bestRsq', z_bestRsq, 'z_mostRF', z_mostRF,  ...
                            'z_mostCh', z_mostCh, 'nanHasMoreCh',nanHasMoreCh , ...
                            'bestIsNotMost', ~ismember(z_bestRsq, z_mostRF) && ~ismember(z_bestRsq, z_mostCh), ...
                            'bestLessChans', max(nChPass(ismember(base_zcos, z_mostCh)) - nChPass(base_zcos == z_bestRsq)));
                        
                        
                        % insert tmpS fields to seriesStruct
                        tmpS_fields = fieldnames(tmpS);
                        for ifield = tmpS_fields'
                            [seriesStruct(:,ibase).(ifield{:})] =  deal(tmpS.(ifield{:}));
                        end
                        clear tmpS tmpS_fields;
                        
                        % insert tmpI fields to expInfo (redundant across zscores)
                        tmpI_fields = fieldnames(tmpI);
                        for ifield = tmpI_fields'
                            seriesStruct(1,ibase).expInfo.(ifield{:}) = tmpI.(ifield{:});
                            [seriesStruct(2:end,ibase).expInfo] = deal(seriesStruct(1,ibase).expInfo);
                        end
                        clear tmpI tmpI_fields;
                    end% ibase
                    
                    
                    % save structure for series
                    serBase = regexp(basenames{serIndices(ibase)}, '(.+)(?=_\()','tokens','once');
                    serBase = serBase{:};
                    save(fullfile(serDir,[serBase '.mat']), 'seriesStruct');
                end % iser
            end %mouse
        end %% makeSeriesStructs
        
        %% makeSumStruct
        
        if makeSumStruct || ~isdir(sumSerDir)
            rng(1); % for reproducibility
            
            if ~isdir(sumSerDir)
                mkdir(sumSerDir);
            end
            
            % get all the series basenames
            tmp = dir(serDir);
            baseseries = {tmp(~cellfun(@isempty, strfind({tmp.name}, '.mat'))).name};
            clear tmp;
            
            % logical of series that dont have any RFs
            badSeries = false(size(baseseries));
            badSeries = struct('baseseries', baseseries, 'badLog', num2cell(badSeries));
            for iser = 1:numel(baseseries)
                % get basic info
                mouse_id = cellstr(regexp(baseseries{iser}, '(.+)(?=_[0-9]+\.mat)', 'tokens', 'once'));
                mouse_id = mouse_id{:};
                series_num = regexp(baseseries{iser}, '([0-9]+)(?=\.mat)', 'tokens', 'once');
                series_num = str2double(series_num{:});
                
                
                tmp = load(fullfile(serDir, baseseries{iser}));
                seriesStruct = tmp.seriesStruct; clear tmp;
                
                % ONLY CLEAN EXPERIMENTS
                if onlyCleanedExps
                    cleanExps = getCleanSparseNoiseExps({mouse_id, series_num});
                    % only need to look at 1 rsq value, bc all the experiments should be the same
                    badExpLog = arrayfun(@(s) any(ismember(cleanExps.badExps, s.exps)), seriesStruct(1,:));
                    seriesStruct = seriesStruct(:,~badExpLog);
                    % sanity check
                    if ~isempty(cleanExps.cleanExps)
                        assert(all(arrayfun(@(s) any(ismember([cleanExps.cleanExps], s.exps)), seriesStruct(:))), ...
                            'Something is wrong with cleanExps. Some experiments are missing');
                    end
                end
                
                % number of experiment types
                nexps = size(seriesStruct,2);
                
                % get the number of channels
                nch = numel(seriesStruct(1).(rfTypeFullString{1}));
                
                % options for fitgmdist
                fitopts = statset('Display', 'off', 'MaxIter', 1e4);
                
                % test to see if any RFs are good
                tmpTest = cellfun(@(rf) [seriesStruct.(rf)], rfTypeFullString, 'uniformoutput', 0);
                tmpTest = cat(2,tmpTest{:});
                if ~any([tmpTest.r2] > rsq_co)
                    badSeries(iser).badLog = true;
                    continue;
                end
                
                % get the zco with the best rsq
                best_zidx = arrayfun(@(e) find([seriesStruct(:,e).bestRsq]), 1:nexps);
                
                % initialize structure
                sumStruct = struct('mouse_id', mouse_id, 'series_num', series_num, 'ch', {[]}, ...
                    'hasRF', false, 'hasBestRF', false,'allGauss', {[]},'bestGauss', {[]}, 'bestLog_all', {[]}, 'nsamp_all',0, 'nsamp_best', 0, ...
                    'gm_all', {[]},'gm_best', {[]}, 'mu_all', {[]}, 'C_all', {[]}, 'mu_best', {[]}, 'C_best',{[]}, 'expInfo_all', {[]});
                
                sumStruct = repmat(sumStruct,1,nch);
                tmpch = num2cell(1:nch);
                [sumStruct.ch] = deal(tmpch{:}); clear tmpch;
                
                S = seriesStruct(:);
                for ich = 1:nch
                    % first test that there is something in the ch
                    chTest = cell2mat(cellfun(@(rf) arrayfun(@(s) s.(rf)(ich).r2 > rsq_co, S), rfTypeFullString, 'uniformoutput',0));
                    if ~any(chTest(:))
                        continue;
                    else
                        sumStruct(ich).hasRF = true;
                    end
                    
                    % get indices into S and rfType for all rfs > rsq_co for the channel
                    [sidx, ridx] = find(chTest);
                    % get the best z idx for channel
                    sidx_bestZ = intersect(find([S.bestRsq]), sidx);
                    
                    % all gaussians with rsq > rsq_co
                    ngauss = numel(sidx);
                    allGauss = arrayfun(@(si, ri) S(si).(rfTypeFullString{ri})(ich), sidx, ridx);
                    % put in x,y coords
                    expInfo_all = ([S(sidx).expInfo]);
                    [expInfo_all.exps] = deal(S(sidx).exps);
                    [allGauss.xcoords] = deal(expInfo_all.xcoords);
                    [allGauss.ycoords] = deal(expInfo_all.ycoords);
                    % and the monitor positions
                    [allGauss.mon_ang] = deal(expInfo_all.mon_ang);
                    [allGauss.mon_elv] = deal(expInfo_all.mon_elv);
                    
                    
                    % put in other fields
                    [allGauss.rfType] = deal(rfTypeFullString{ridx});
                    [allGauss.exps] = deal(S(sidx).exps);
                    [allGauss.z_cutoff] = S(sidx).z_cutoff;
                    
                    % get new pars and cov matrices that are scaled to one
                    [newImg, newC] = arrayfun(@(g) TwoDGaussian(g.pars,{g.xcoords, g.ycoords}, true, any(strfind(thisFit, 'chol'))),  allGauss, 'uniformoutput', 0);
                    [allGauss.fitImg] = deal(newImg{:});
                    [allGauss.C] = deal(newC{:});
                    
                    
                    % see if there was a bad monitor position if multiple
                    % experiments
                    
                    chExps = unique(cellfun(@(i) num2str(i), {allGauss.exps}, 'uniformoutput',0));
                    % need to do it this way via strings for some reason?
                    tmp  = cellfun(@strsplit, chExps, 'uniformoutput', 0);
                    chExps = cellfun(@str2double, tmp, 'uniformoutput',0);
                    % exps where the ch had RFs
                    nexps4Ch = numel(chExps);
                    
                    
                    
                    if nexps4Ch  > 1
                        allMonPos = [[allGauss.mon_ang];[allGauss.mon_elv]]';
                        unqMonPos = unique(allMonPos, 'rows');
                        nMonPos = size(unqMonPos,1);
                        
                        allStimPos = [[expInfo_all.stim_xc]; [expInfo_all.stim_yc]]';
                        unqStimPos = unique(allStimPos, 'rows');
                        nStimPos = size(unqStimPos, 1);
                        % if the experiments change bc monitor positions is changed
                        if nMonPos == nexps4Ch
                            useMonitorPos = true;
                        elseif  nStimPos == nexps4Ch
                            useMonitorPos = true;
                            allMonPos = allStimPos;
                            unqMonPos = unqStimPos;
                            nMonPos = nStimPos;
                        else
                            useMonitorPos = true;
                        end
                        
                        all_xy = [[allGauss.x0]; [allGauss.y0]]';
                        
                        % multiple monitor positions
                        if useMonitorPos
                            % pairwise distances
                            pjitter = squareform(pdist(all_xy));
                            % see the jitter btw pairs within a monitor position
                            wInPos = arrayfun(@(ipos) squareform(pjitter(ismember(allMonPos, unqMonPos(ipos,:), 'rows'),ismember(allMonPos, unqMonPos(ipos,:), 'rows')), 'tovector'), 1:nMonPos, 'uniformoutput',0);
                            %wInMed = cellfun(@median, wInPos);
                            wInPos = cellfun(@(i) i(:)', wInPos, 'uniformoutput', 0);
                            wInMed = median([wInPos{:}]);
                            
                            % across positions
                            axPos = arrayfun(@(ipos) unique(pjitter(ismember(allMonPos, unqMonPos(ipos,:), 'rows'),~ismember(allMonPos, unqMonPos(ipos,:), 'rows'))), 1:nMonPos, 'uniformoutput', 0);
                            axPos = cellfun(@(i) i(:)', axPos, 'uniformoutput', 0);
                            %axMed = cellfun(@median, axPos);
                            axMed = median([axPos{:}]);
                            
                            if axMed > 1.5*wInMed || axMed > 8 || wInMed > 8
                                fprintf('\t check iser%d,ch%d; axMed = %0.1f, wInMed = %0.1f \n', iser, ich, axMed , wInMed);
                                if wInMed > 8
                                    warning('check iser%d,ch%d; axMed = %0.1f, wInMed = %0.1f \n', iser, ich, axMed , wInMed);
                                end
                                if axMed > 8
                                    fprintf('\n \n rsq = %d, iser = %d, ch = %d, axMed = %0.1f, wInMed = %0.1f, \n \n', rsq_co * 100, iser, ich, axMed,wInMed )
                                end
                            end
                            
                        else
                            warning('nexp ~= nMonPos/nStimPos! iser = %d', iser);
                        end
                    end
                    
                    
                    % gaussians with the best z_co
                    bestLog_all = ismember(sidx, sidx_bestZ);
                    bestGauss = allGauss(bestLog_all);
                    %%%%%%%%%%%%%%%%%%%%%%
                    %%%% FOR allGauss %%%%
                    %%%%%%%%%%%%%%%%%%%%%%
                    % created weighted samples
                    w_all = [allGauss.r2]/sum([allGauss.r2]);
                    % number of samples to generate for each gauss
                    nsamp_tot = 1e5;
                    nsamp_all = round(w_all*nsamp_tot);
                    % produce samples by doing a cholesky decomp (there is only 1 of these per matrix)
                    % L = chol(C);
                    % samples (x) : x = Ly + u where y ~ N(0,1);
                    all_samples = arrayfun(@(g,smp) chol(g.C, 'lower')*normrnd(0,1, 2,smp) + repmat([g.x0; g.y0], 1, smp), allGauss, nsamp_all', 'uniformoutput', 0);
                    all_samples = cell2mat(all_samples');
                    % fit new gaussian; k = 1 (want 1 gaussian back)
                    gm_all = fitgmdist(all_samples',1, 'RegularizationValue',0.01,'CovarianceType','full', ...
                        'Replicates', 10, 'SharedCovariance', false, 'options', fitopts);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%
                    %%%% FOR bestGauss %%%%
                    %%%%%%%%%%%%%%%%%%%%%%%
                    if any(bestLog_all)
                        % created weighted samples
                        w_best = [bestGauss.r2]/sum([bestGauss.r2]);
                        
                        % number of samples to generate for each gauss
                        nsamp_best = round(w_best*nsamp_tot);
                        % produce samples
                        best_samples = arrayfun(@(g,smp) chol(g.C, 'lower')*normrnd(0,1, 2,smp) + repmat([g.x0; g.y0], 1, smp), bestGauss, nsamp_best', 'uniformoutput', 0);
                        best_samples = cell2mat(best_samples');
                        % fit new gaussian
                        gm_best = fitgmdist(best_samples',1, 'RegularizationValue',0.01,'CovarianceType','full', ...
                            'Replicates', 10, 'SharedCovariance', false, 'options', fitopts);
                        % assign to structure
                        sumStruct(ich).hasBestRF = true;
                        sumStruct(ich).nsamp_best = nsamp_best;
                        sumStruct(ich).gm_best = gm_best;
                        sumStruct(ich).mu_best = gm_best.mu;
                        sumStruct(ich).C_best = gm_best.Sigma;
                    end
                    
                    % assign to structure
                    sumStruct(ich).allGauss = allGauss;
                    sumStruct(ich).bestGauss = bestGauss;
                    sumStruct(ich).bestLog_all = bestLog_all;
                    sumStruct(ich).nsamp_all = nsamp_all;
                    sumStruct(ich).gm_all = gm_all;
                    sumStruct(ich).mu_all = gm_all.mu;
                    sumStruct(ich).C_all = gm_all.Sigma;
                    sumStruct(ich).expInfo_all = expInfo_all;
                    
                end % ich
                
                % if v1 average together all the channels and get 1 gauss
                % for all channels
                if runV1
                    % initialize
                    avgV1Rf = struct('mouse_id', mouse_id, 'series_num', series_num, 'all_chans', {[]}, ...
                        'best_chans', {[]}, 'hasRF', true, 'hasBestRF', false, 'nsamp_all',0, 'nsamp_best', 0, ...
                        'gm_all', {[]},'gm_best', {[]}, 'mu_all', {[]}, 'C_all', {[]}, 'mu_best', {[]}, 'C_best',{[]});
                    
                    allCh = sumStruct([sumStruct.hasRF])';
                    bestCh = sumStruct([sumStruct.hasBestRF])';
                    if onlyCleanedExps
                        allCh = allCh(~ismember([allCh.ch],[cleanExps.badChans]));
                        bestCh = bestCh(~ismember([bestCh.ch],[cleanExps.badChans]));
                    end
                    
                    % AVG ACROSS CHANNELS
                    % created weighted samples
                    w_all = ones(1, numel(allCh)) ./ numel(allCh);
                    % number of samples to generate for each gauss
                    nsamp_tot = 1e5;
                    nsamp_all = round(w_all*nsamp_tot);
                    % produce samples by doing a cholesky decomp (there is only 1 of these per matrix)
                    % L = chol(C);
                    % samples (x) : x = Ly + u where y ~ N(0,1);
                    all_samples = arrayfun(@(g,smp) chol(g.C_all, 'lower')*normrnd(0,1, 2,smp) + repmat(g.mu_all', 1, smp), allCh, nsamp_all', 'uniformoutput', 0);
                    all_samples = cell2mat(all_samples');
                    % fit new gaussian; k = 1 (want 1 gaussian back)
                    gm_all = fitgmdist(all_samples',1, 'RegularizationValue',0.01,'CovarianceType','full', ...
                        'Replicates', 10, 'SharedCovariance', false, 'options', fitopts);
                    
                    % AVG ACROSS CHANNELS - BEST Z-Score
                    if numel(bestCh)
                        % created weighted samples
                        w_best = ones(1, numel(bestCh)) ./ numel(bestCh);
                        % number of samples to generate for each gauss
                        nsamp_tot = 1e5;
                        nsamp_best = round(w_best*nsamp_tot);
                        % produce samples by doing a cholesky decomp (there is only 1 of these per matrix)
                        % L = chol(C);
                        % samples (x) : x = Ly + u where y ~ N(0,1);
                        best_samples = arrayfun(@(g,smp) chol(g.C_best, 'lower')*normrnd(0,1, 2,smp) + repmat(g.mu_best', 1, smp), bestCh, nsamp_best', 'uniformoutput', 0);
                        best_samples = cell2mat(best_samples');
                        % fit new gaussian; k = 1 (want 1 gaussian back)
                        gm_best = fitgmdist(best_samples',1, 'RegularizationValue',0.01,'CovarianceType','full', ...
                            'Replicates', 10, 'SharedCovariance', false, 'options', fitopts);
                        
                        % assign to structure
                        avgV1Rf.hasBestRF = true;
                        avgV1Rf.nsamp_best = nsamp_best;
                        avgV1Rf.gm_best = gm_best;
                        avgV1Rf.mu_best = gm_best.mu;
                        avgV1Rf.C_best = gm_best.Sigma;
                        avgV1Rf.best_chans =[bestCh.ch];
                    end % numel(bestCh)
                    
                    % assign to structure
                    avgV1Rf.nsamp_all = nsamp_all;
                    avgV1Rf.gm_all = gm_all;
                    avgV1Rf.mu_all = gm_all.mu;
                    avgV1Rf.C_all = gm_all.Sigma;
                    avgV1Rf.all_chans =[allCh.ch];
                    
                    % save new structure;
                    [~,fname, ~] = fileparts(baseseries{iser});
                    fname = fullfile(sumSerDir,[fname '_avgV1Rf.mat']);
                    save(fname, 'avgV1Rf'); clear fname
                end % runv1 (gets "avg ch" from all the ch averages)
                
                % save new structure;
                fname = fullfile(sumSerDir,[baseseries{iser}]);
                save(fname, 'sumStruct'); clear fname;
            end % iser
            
            % save badLog
            if any([badSeries.badLog])
                fname = fullfile(sumSerDir, 'badSeries.mat');
                save(fname, 'badSeries');clear fname;
            end
        end %% makeSumStruct
    end % for rsq
end % if makeSumStruct || makeSeriesStruct



%% plotSeriesStruct
if plotSeriesStruct
    setFigShit
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% PLOT SERIES STRUCT %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compares rsq_co?
    RF_directoryFunctions; % run script (same as what is commented out below)
%     % inline functions for easier access to rsq_co folder
%     % rfCdir : directory for summed gaussian structs (depends on rsq)
%     % sumDir : directory that contains sum structs and series structs folders
%     % sumSerDir : directory for sum structs
%     % serDir : directory for series structs
%     get_rfCdir = @(my_rsq) fullfile(muaDir, sprintf('summedFits_rsq%d', my_rsq*100));
%     get_sumDir = @(my_rsq) fullfile(get_rfCdir(my_rsq), thisFit);
%     if onlyCleanedExps
%         get_sumSerDir = @(my_rsq) fullfile(get_sumDir(my_rsq),'cleanExpFits');
%     else
%         get_sumSerDir = @(my_rsq) fullfile(get_sumDir(my_rsq),'seriesFits');
%     end
%     get_serDir = @(my_rsq) fullfile(get_sumDir(my_rsq), 'seriesInfo');
%     
%     % same, but for z_co folders for plotting
%     % zDir : has RFInfo, for plotting the interpolated data
%     %       e.g. RF = load(fullfile(zDir, [basenames{ibase},'_RFInfo.mat']));
%     get_zDir = @(my_zco) fullfile(muaDir, sprintf('z_co%d',my_zco));
    
    
    % get all the series basenames (not all will be present in sum structs but
    % will have a series struct regardless of rsq_co
    tmp = dir(get_serDir(0.45));
    baseseries = {tmp(~cellfun(@isempty, strfind({tmp.name}, '.mat'))).name};
    clear tmp;
    
    % plot different cutoffs
    plot_rsqs = [0.45 0.5 0.6]; %here!
    if any(setdiff(myRsqCutoffs,plot_rsqs))
        warning('not going to plot rsq %s. change me here.', strjoin(arrayfun(@num2str,setdiff(myRsqCutoffs,plot_rsqs), 'uniformoutput',0),','));
    end
    %rsq_colors = {'c', 'g', 'r'};
    rsq_colors = {blue, green, red};
    % difft linewidhts just to see better
    rsq_lw = (2:-0.25:1.5);
    % v1 avg color
    v1_cols = {purple,orange,ochre};
    v1_markers = {'o', 's'};
    % adjust above based off myRsqCutoffs
    tmpidx = ismember(plot_rsqs,myRsqCutoffs);
    if ~all(tmpidx) 
        plot_rsqs = plot_rsqs(tmpidx); rsq_colors = rsq_colors(tmpidx); rsq_lw = rsq_lw(tmpidx); v1_cols = v1_cols(tmpidx);
    end; clear tmpidx
    
    % inline function for shading
    shadeColor = @(col,idx) (nonzeros(mapMyColor(col) .* sparse(idx*ones(1,3), 1:3, ones(1,3), 200,3))');
    
    
    % gaussian elipse alpha
    my_alpha = 0.75;
    
    for iser = 1:numel(baseseries)
        
        % get basic info
        mouse_id = cellstr(regexp(baseseries{iser}, '(.+)(?=_[0-9]+\.mat)', 'tokens', 'once'));
        mouse_id = mouse_id{:};
        series_num = regexp(baseseries{iser}, '([0-9]+)(?=\.mat)', 'tokens', 'once');
        series_num = str2double(series_num{:});
        
        %index in basenames
        baseIndices = find(strcmp(allMice, mouse_id) & allSeries == series_num);
        % experiments in basenames (note: all may not be used)
        [~,baseExps] = arrayfun(@(i) ....
            regexp(basenames{i}, '\(([0-9,\,]+)\)', 'match', 'tokens', 'once'),...
            baseIndices, 'uniformoutput', 0);
        baseExps = cellfun(@(ibase) cellfun(@str2num,strsplit(ibase, ',')), ...
            [baseExps{:}], 'uniformoutput',0);
        
        if onlyCleanedExps
            cleanExps = getCleanSparseNoiseExps({mouse_id, series_num});
            badExpLog = cellfun(@(ibase) any(ismember(ibase, cleanExps.badExps)),baseExps);
            baseExps = baseExps(~badExpLog);
            baseIndices = baseIndices(~badExpLog);
            % sanity check
            if ~isempty(cleanExps.cleanExps)
                assert(all(cellfun(@(base) any(ismember([cleanExps.cleanExps], base)), baseExps)), ...
                    'Something is wrong with cleanExps. Some experiments are missing');
            end
        end
        nbaseExps = numel(baseExps);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%% LOAD SUM STRUCTS %%%
        %%%%%%%%%%%%%%%%%%%%%%%%
        % for different rsqs... check that files exist
        rsqLog = logical(arrayfun(@(rsq) exist(fullfile(get_sumSerDir(rsq), baseseries{iser}), 'file'), plot_rsqs));
        
        if ~any(rsqLog)
            warning('Can not locate any files for %s', baseseries{iser});
            continue
        end
        
        % rsqs
        my_rsqs = plot_rsqs(rsqLog);
        my_cols = rsq_colors(rsqLog);
        my_lw = rsq_lw(rsqLog);
        nrsqs = nnz(rsqLog);
         
        % load all the files
        S = arrayfun(@(rsq) load(fullfile(get_sumSerDir(rsq), baseseries{iser})), my_rsqs);
        % number of channels
        nch = numel(S(1).sumStruct);
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% V1 only: get RF across chans %%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if runV1
            [~,fname, ~] = fileparts(baseseries{iser});
            fname = [fname '_avgV1Rf.mat'];  %#ok<AGROW>
            V1 = arrayfun(@(rsq) load(fullfile(get_sumSerDir(rsq), fname)), my_rsqs); 
            clear fname
            tmp = num2cell(my_rsqs);
            [V1.rsq] = deal(tmp{:}); clear tmp
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% Z IS BEST ? %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % by best z_cutoff or across z_cutoffs?
        for zIsBest = [false true]
            
            if zIsBest
                mu = 'mu_best';
                C = 'C_best';
                gauss = 'bestGauss'; % indv fits
                zTitle = 'zBEST_ONLY';
                v1chans = 'best_chans'; % for avg v1
            else
                mu = 'mu_all';
                C = 'C_all';
                gauss = 'allGauss'; % indv fits
                zTitle = 'ALL_Z';
                v1chans = 'all_chans'; % for avg v1
            end
            % logs of channels present
            % [rsq x ch]
            if zIsBest
                rsqChLogs = cell2mat(arrayfun(@(s) [s.sumStruct.hasBestRF], S, 'uniformoutput',0)');
            else
                rsqChLogs = cell2mat(arrayfun(@(s) [s.sumStruct.hasRF], S, 'uniformoutput',0)');
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%
            % inline function to index S by rsq
            % missing channels will not show up
            if exist('getSbyRsq', 'var'); clear rsqS; end;
            getSbyRsq = @(rsq) S(my_rsqs == rsq).sumStruct(rsqChLogs(my_rsqs== rsq,:));
            %%%%%%%%%%%%%%%%%%%%%
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%% KEEP TRACK OF WHAT EXPERIMENTS WERE USED?  %%%%%%%%%%%%%
            %%%%%% this is harder for the allGauss case bc    %%%%%%%%%%%%%
            %%%%%% each channel may be different based on fit %%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%% need it for  best Z for the exp/rsq_co      %%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%% zforRsq has nExps x nRsq %%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            zForRsq = struct('rsq', repmat(num2cell(my_rsqs), nbaseExps,1), ...
                'exps', cell(nbaseExps,nrsqs), 'z_best',cell(nbaseExps,nrsqs));
            
            for irsq = 1:nrsqs
                if zIsBest
                    % uses inline function getSbyRsq
                    tmp = cell2mat(arrayfun(@(s) s.expInfo_all, ...
                        getSbyRsq(my_rsqs(irsq)), 'uniformoutput',0));
                    % tmpidx is the index to unique instances of experiment
                    % (only 1 best z_co)
                    [~,tmpidx, ~] = unique(cellfun(@(x)(mat2str(x)),{tmp.exps},'uniformoutput',false));
                    zrsq = {tmp(tmpidx).z_bestRsq};
                else
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%% at least see if any z_cos/exp not used at all? %%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    tmp = cell2mat(arrayfun(@(s) s.allGauss, getSbyRsq(my_rsqs(irsq)), 'uniformoutput',0)');
                    % tmpass is where each unq exp was: across all relevant channels
                    % ie chans with RFs that passed the rsq_co
                    % (multiple zco for exp)
                    [~,tmpidx, tmpass] = unique(cellfun(@(x)(mat2str(x)),{tmp.exps},'uniformoutput',false));
                    % get which zcos were considerd
                    tmpu = unique([tmpass, cat(1, tmp.z_cutoff)], 'rows');
                    zrsq = arrayfun(@(iexp) tmpu(tmpu(:,1) == iexp,2)', 1:numel(tmpidx), 'uniformoutput',0);
                end % zIsBest
                zexps = {tmp(tmpidx).exps};
                nrsqExps = numel(zrsq);
                zForRsq(1:nrsqExps,irsq) = cell2struct([repmat({my_rsqs(irsq)},1,nbaseExps); zexps; zrsq], {'rsq', 'exps', 'z_best'});
                clear tmp tmpidx tmpass
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%% PLOT EACH EXPERIMENT AND SEE HOW WELL FITS %%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for ibase = 4:nbaseExps
                
                % base name
                ibasename = basenames{baseIndices(ibase)};
                % which z_cos are involved?
                zlog = any((arrayfun(@(z) any(ismember(z.exps , baseExps{ibase})), zForRsq)),2);
                if ~any(zlog)
                    warning('no good RFs for %s', ibasename);
                end
                
                % get the relevant z_cutoffs
                base_zcos = unique([zForRsq(zlog,:).z_best]);
                
                % load all the RFs for base_zcos
                zRFs = arrayfun(@(zco) load(fullfile(get_zDir(zco),[ibasename,'_RFInfo.mat'])), base_zcos);
                
                for izco = 1:numel(base_zcos)
                    RF = zRFs(izco);
                    for irf = 1:nrf
                        thisRF = rfTypeFullString{irf};
                        
                        figure(figPars, 'Name', sprintf('%s,z_co=%d %s %s',zTitle, base_zcos(izco), ibasename, thisRF), ...
                            'position', [rfPos(irf) zPos(izco), scrnfigw, scrnfigh])
                        
                        h = plotRFbyChStruct(RF.(thisRF), RF.sRegionX, RF.sRegionY,RF.intRegionX,RF.intRegionY,subplotsAreTight);
                        
                        for ich = 1 : nch
                            
                            for irsq = 1:nrsqs
                                % DATA
                                tmpS = getSbyRsq(my_rsqs(irsq));
                                chLog = [tmpS.ch] == ich;
                                if ~any(chLog)
                                    continue
                                else
                                    hold(h(ich), 'on')
                                end
                                
                                % did the channel in this exp get used during
                                % the averaging? (will only matter for the individual fits)
                                chInExp = arrayfun(@(s) any(ismember(s.exps, baseExps{ibase})), tmpS(chLog).(gauss))';
                                if isempty(chInExp)
                                    if zIsBest && ~all([tmpS(chLog).hasBestRF])
                                        continue
                                    else
                                        error('whyyy')
                                    end
                                end
                                
                                % INDIVIDUAL FITS
                                indvLog = chInExp & [tmpS(chLog).(gauss).z_cutoff] == base_zcos(izco) & ...
                                    strcmp({tmpS(chLog).(gauss).rfType}, thisRF);
                                
                                if any(indvLog)
                                    indvFit = tmpS(chLog).(gauss)(indvLog);
                                    indvCol = shadeColor(my_cols{irsq},180);
                                    ilw = 0.75;
                                    
                                    plot(h(ich), indvFit.x0,indvFit.y0,...
                                        'color', indvCol , 'marker', '+', 'linestyle', 'none', 'linewidth',ilw);
                                    
                                    plot(h(ich), indvFit.x_fit, indvFit.y_fit,  'color', indvCol ,...
                                        'marker', 'none', 'linestyle', '-','linewidth', ilw);
                                end
                                
                                
                                % OVERALL AVERAGE FIT
                                [x,y] = getGaussEllipse(tmpS(chLog).(mu)', tmpS(chLog).(C), my_alpha);
                                plot(h(ich), tmpS(chLog).(mu)(1),tmpS(chLog).(mu)(2),...
                                    'color', my_cols{irsq} , 'marker', '+', 'linestyle', 'none', 'linewidth', my_lw(irsq));
                                plot(h(ich), x,y,  'color', my_cols{irsq} , 'marker', 'none', 'linestyle', '-','linewidth', my_lw(irsq));
                                
                                
                                
                                 
                                % V1 AVG CHANNEL
                                if runV1
                                    v1Log = [V1.rsq] == my_rsqs(irsq);
                                    thismarker = v1_markers{any(V1(v1Log).avgV1Rf.(v1chans)==ich) + 1};
                                    [x,y] = getGaussEllipse(V1(v1Log).avgV1Rf.(mu)', V1(v1Log).avgV1Rf.(C), my_alpha);
                                    plot(h(ich), V1(v1Log).avgV1Rf.(mu)(1),V1(v1Log).avgV1Rf.(mu)(2),...
                                        'color', v1_cols{irsq} , 'marker', thismarker, 'linestyle', 'none', 'linewidth', my_lw(irsq));
                                    plot(h(ich), x,y,  'color', v1_cols{irsq} , 'marker', 'none', 'linestyle', '-','linewidth', my_lw(irsq));
                                end
                                
                            end % nrsqs
                        end %ich
                    end % irf
                end % izco
                pause
            end % ibase
        end % zIsBest
        pause; close all;
    end % iser
    
end % plotSeriesStruct