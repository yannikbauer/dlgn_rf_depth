function sparseNoiseExpByChan(key, RFtwin, zcutoff, file_ext, plotFlag, figDir, matDir, saveRaw) %, loadMatFile)
% defulat saveRaw = false, file_ext can be maue or filtdat (2K filtered dat)
% todo: should have option for reading matfile (either rawEvents or RFinfo)
% todo: plotting position on screen!!! 

if nargin < 2 || isempty(RFtwin)
    RFtwin = [0 0.35];
end

if nargin < 3 || isempty(zcutoff)
    zcutoff = 8;
end

if nargin < 4 || isempty(file_ext)
    % file_ext = 'muae';
    file_ext = 'envl.dat';
end

if ~strcmpi(file_ext,'muae') && ~strcmpi(file_ext,'envl.dat')
    error('other than reading your %s files as an LFP, not all conds implemented')
end

if nargin < 5 || isempty(plotFlag)
    plotFlag = false;
end

if nargin < 6  || isempty(figDir)
    figDir = [];
else
    figDir = fullfile(figDir, sprintf('z_cutoff_%d/', zcutoff)); 
    if ~isdir(figDir); mkdir(figDir); end
end

if nargin < 7  || isempty(matDir)
    matDir = [];
else
    matDir = fullfile(matDir, sprintf('z_cutoff_%d/', zcutoff)); 
    if ~isdir(matDir);mkdir(matDir);end
end

if nargin < 8 || isempty(saveRaw)
    saveRaw = false;
elseif saveRaw && isempty(matDir)
    error('need matDir for saveRaw');
else
    rawDir = fileparts(matDir);
end

% fetch mouse id if you must....
if plotFlag || ~isempty(figDir) || ~isempty(matDir)
    mouse_id = fetch1(data.Mice(key(1)), 'mouse_id');
end


% define variables
% 2d Gaussian filter
sigma = 1;

% for smoothing/filtering RFs after interpolation (ASK LAURA)
cutoff = ceil(3*sigma);
h = fspecial('gaussian',2*cutoff+1,sigma);

% experiments as string
expStr = ['(' sprintf('%d,', [key.exp_num])];
expStr(end) = ')';

% number of exps
nexp = numel(key);

% get array layout and the number of channels
% array layout
array_layout = fetch1(data.Arrays & data.Series(key(1)), 'array_layout');

nch = numel(array_layout);

% read lfp or muae and keep track of experiment times
% note that for MUAE : the channels are read from the nsx file so the array
% layout matters
% does that mean for lfp, the chans are diff? (ASK LAURA)
myData = [];
t_myData = [];
eofTime = nan(1, nexp);

for ikey = 1 : nexp
    try
        if strcmpi(file_ext,'muae')
            [temp_lfp, temp_tlfp] = readLFP(key(ikey).mouse_counter, key(ikey).series_num, key(ikey).exp_num, 'data', file_ext);
            sample_rate = 1250;
        elseif strcmpi(file_ext,'envl.dat')
            % get path to nsx files
            [unsortedFolder, ~, sortedDataFolder] = getPathTo(key(ikey), 'data');
            
            % directory where envl.dat files are
%             datDir =  '~/busseDB/Sinem/'; % TK: will need to be changed
                                          % datDir = sortedDataFolder; 
            datDir =  '/Volumes/lab/users/yannik/projects/ret2dlgn/analyses/dlgn_rf_depth/data/';%'/Volumes/lab/users/miro/data/RF/spyke'; % TK: will need to be changed

            dtemp = dir(datDir);
            % regexp pattern for .dat file
            pattern = [mouse_id,'_', sprintf('s%02d',key(ikey).series_num),'.+', ...
                sprintf('%02d.', key(ikey).exp_num), 'ns6.', file_ext, '$'];
            % get the file name
            fnameLog = ~cellfun(@isempty, regexp({dtemp.name},pattern, 'match'));
            assert(nnz(fnameLog)==1, 'cant find file or too many files found. \n \t using regexp pattern %s', pattern);
            fname = dtemp(fnameLog).name;
            myDat = fullfile(datDir,fname);
            % read json file
            myJson = [myDat,'.json'];
            assert(exist(myJson, 'file') > 0,['cant find ' myJson]);
            fid = fopen(myJson);
            raw = fread(fid,inf);
            str = char(raw');
            fclose(fid);
            info = JSON.parse(str);
            clear fid raw str
            % make sure the offsets match up
            [eofTimes_nsx, startDelay] = nsxEOF(unsortedFolder);
            offsetsmatch= ceil(startDelay{key(ikey).exp_num}/30000 * info.sample_rate) == info.nsamples_offset;
            assert(offsetsmatch, 'offsets between json and nsx file do not match');
            % get the number of samples
            nsamples = floor((eofTimes_nsx(key(ikey).exp_num) - startDelay{key(ikey).exp_num})/ (30000/info.sample_rate));
            % read .dat file
            fid = fopen(myDat,'r');
            temp_lfp = fread(fid,[info.nchans,nsamples], 'int16');
            fclose(fid); clear fid;
            % time and sample rate
            temp_tlfp = ((0:1:nsamples) + info.nsamples_offset)/info.sample_rate;
            sample_rate = info.sample_rate;
        end
    catch ME
        fprintf('Error reading %s for mouse %s, series %02d, exp %02d\n', file_ext, mouse_id, key(ikey).series_num, key(ikey).exp_num);
        fprintf('%s', getReport(ME));
        fprintf('\n');
        RFcoord = []; %#ok<NASGU>
        valRF = []; %#ok<NASGU>
        gradDeg = []; %#ok<NASGU>
        return
    end
    temp_lfp = squeeze(TensorBPTime(3, inf, sample_rate, reshape(temp_lfp, [1 size(temp_lfp)])));
    myData = [myData temp_lfp]; %#ok<AGROW>
    if ikey == 1
        t_myData = temp_tlfp;
    else
        t_myData = [t_myData temp_tlfp + t_myData(end)]; %#ok<AGROW>
    end
    eofTime(ikey) = t_myData(end);
end


% get the stimulus onsets
t_info = fetch(data.GratingTrials(key), 'grat_num', 'trial_onset', 'trial_offset');
if isempty(t_info)
    error('singleExpSparseNoise:noGratingTrials', 'Populate data.Stimuli for mc %d s%02d e%02d', ...
        key(1).mouse_counter, key(1).series_num, [key.exp_num]);
end
[t_info.trial_onset] = vecdeal([t_info.trial_onset] / 30000);
[t_info.trial_offset] = vecdeal([t_info.trial_offset] / 30000);


% sort trials by exp and trial_num
if ~issorted([t_info.trial_num])
    [~, sidx] = sortrows([[t_info.exp_num];[t_info.trial_num]]');
    t_info = t_info(sidx);
end

% adjust the stimulus onsets and offsets for concatenated experiments
for ikey = 2 : nexp
    expLog = [t_info.exp_num] == key(ikey).exp_num;
    [t_info(expLog).trial_onset]  = vecdeal([t_info(expLog).trial_onset]  + eofTime(ikey-1));
    [t_info(expLog).trial_offset] = vecdeal([t_info(expLog).trial_offset] + eofTime(ikey-1));
end

% grating numbers and total number of gratings
grats = unique([t_info.grat_num]);
ngrats = numel(grats);

% get the spatial configuration of gratings
% we can use key(1) because if there are more keys then they would refer to
% identical experiments
[g_info, active_par_names] = fetchn(data.StimInfo(key(1)), ...
    'num_active_grats', 'active_par_names');
g_info = g_info{1};

active_par_names = active_par_names{1};
% check the order of parameters
if ~strcmp(active_par_names{1}, 'grat_y_position') || ...
    ~strcmp(active_par_names{2}, 'grat_x_position') || ...
    ~strcmp(active_par_names{3}, 'grat_contrast')
    error('sparseNoiseLFP:oldFormat', 'Cannot deal with the way stimInfo is organized. Repopulate data.StimInfo for %s, s%02d, e%02d', ...
        mouse_id, key(1).series_num, key(1).exp_num)
end



% figure out on/off gratings and their idx in g_info
if fetch1(data.GratingConditions(key(1), sprintf('grat_num = %d', g_info(1,1,1))), 'grat_contrast')==1
    onStimIdx = 1; offStimIdx = 2;
else
    onStimIdx = 2; offStimIdx = 1;
end


% get the individual events for each trial

% number of reps per grat
nreps = numel(t_info)/ngrats;
% preallocate IndividualEvents
IndividualEvents = cell(1, ngrats);
eventCh = cell(1, ngrats);
eventTidx = cell(1, ngrats);

fprintf('\n \t calculating individual events ... \t ')
tic
for igrat = 1 : ngrats
    trial_log = [t_info.grat_num] == igrat;
    [~, et, temp] = TensorEvent(reshape(myData,...
        [1 size(myData)]), [t_info(trial_log).trial_onset], [-0.1 0.4], t_myData(end), 'meansub');
    % channel x time
    IndividualEvents{igrat} = squeeze(temp);
    temp1 = arrayfun(@(ich) squeeze(IndividualEvents{igrat}(ich,:,:))', 1:nch, 'uniformout', 0);
    IndividualEvents{igrat} = cat(1,temp1{:});
    eventCh{igrat} = repelem((1:nch)', nreps,1);
    eventTidx {igrat} = repmat(find(trial_log(:)),nch,1);
end
toc

allEvents = cat(1,IndividualEvents{:});
eventCh = cat(1,eventCh{:});
eventTidx = cat(1,eventTidx{:});

% z-score and find elements larger than cutoff
zevents = (allEvents - mean(allEvents(:))) / std(allEvents(:),1);
badz_log = any(abs(zevents) > zcutoff,2);

% also zscore by each channel individually 
zevents_zbych = nan(size(zevents));
for ich = 1:nch
    zevents_zbych(eventCh == ich, :) = ...
        (allEvents(eventCh == ich,:) - mean(reshape(allEvents(eventCh == ich,:),1,[])))/std(reshape(allEvents(eventCh == ich,:),1,[]),1);
end
badz_log_zbych = any(abs(zevents_zbych) > zcutoff,2);


% save raw data as well as badz_log to indicate what was thrown out
if saveRaw        
    fname = fullfile(rawDir , sprintf('%s_%02d_%s_rawEvents.mat', mouse_id, key(1).series_num, expStr));
    fprintf('\n \t saving raw mat %s... \t', fname)
    tic
    % serialize data before daving
    % see http://undocumentedmatlab.com/blog/serializing-deserializing-matlab-data
    allEvents_byteStream = getByteStreamFromArray(allEvents);%#ok<NASGU>
    t_info_byteStream = getByteStreamFromArray(t_info);%#ok<NASGU>
    eventTidx_byteStream = getByteStreamFromArray(eventTidx);%#ok<NASGU>
    eventCh_byteStream = getByteStreamFromArray(eventCh);%#ok<NASGU>
    badz_log_byteStream = getByteStreamFromArray(badz_log);%#ok<NASGU>
    % USE getArrayFromByteStream for deserializing 
    save(fname, ...
        'onStimIdx', 'offStimIdx', 'g_info', 't_info_byteStream', ...
        'allEvents_byteStream', 'eventCh_byteStream', 'eventTidx_byteStream', 'nch', 'badz_log_byteStream','et', '-v6');
    toc
    % return
end




% zscore again without badz trials
zevents = allEvents(~badz_log,:);
zmean = mean(zevents(:));
zstd = std(zevents(:),1);

% note that cleanEvents can be changed so that it has bad trials but is scaled by clean data
cleanEvents = (allEvents(~badz_log,:) - zmean) /zstd;
cleanCh = eventCh(~badz_log);
cleanTidx = eventTidx(~badz_log);


% zscore bych again without badz trials
zevents_zbych = allEvents(~badz_log_zbych,:);
cleanCh_zbych = eventCh(~badz_log_zbych);
cleanTidx_zbych = eventTidx(~badz_log_zbych);
for ich = 1:nch
    zevents_zbych(cleanCh_zbych == ich, :) = ...
        (zevents_zbych(cleanCh_zbych == ich,:) - ...
        mean(reshape(zevents_zbych(cleanCh_zbych == ich,:),1,[])))/...
        std(reshape(zevents_zbych(cleanCh_zbych == ich,:),1,[]),1);
end





% avg across trials and divide into tiles
ny = size(g_info,1); %xpos
nx = size(g_info,2); %ypos
onRespAxTime = cell(nx, ny);
offRespAxTime = cell(nx, ny);
onRespAxTime_zbych = cell(nx, ny);
offRespAxTime_zbych = cell(nx, ny);
fprintf('\n \t averaging across trials for on/off stims and tiling... \t ')
tic
for ir = 1:nx
    for ic = 1:ny
        % for regular zscoring 
        onLog = ismember(cleanTidx, find([t_info.grat_num] == g_info(ir,ic,onStimIdx)));
        offLog = ismember(cleanTidx, find([t_info.grat_num] == g_info(ir,ic,offStimIdx)));
        onRespAxTime{ir, ic} = cell2mat(arrayfun(@(ich) nanmean(cleanEvents(onLog & cleanCh == ich, :),1), 1:nch, 'uniformoutput',0)');
        offRespAxTime{ir, ic} = cell2mat(arrayfun(@(ich) nanmean(cleanEvents(offLog & cleanCh == ich, :),1), 1:nch, 'uniformoutput',0)');
        
        % for zscoring by channel
        onLog_zbych = ismember(cleanTidx_zbych, find([t_info.grat_num] == g_info(ir,ic,onStimIdx)));
        offLog_zbych = ismember(cleanTidx_zbych, find([t_info.grat_num] == g_info(ir,ic,offStimIdx)));
        onRespAxTime_zbych{ir, ic} = cell2mat(arrayfun(@(ich) nanmean(zevents_zbych(onLog_zbych  & cleanCh_zbych == ich, :),1), 1:nch, 'uniformoutput',0)');
        offRespAxTime_zbych{ir, ic} = cell2mat(arrayfun(@(ich) nanmean(zevents_zbych(offLog_zbych & cleanCh_zbych == ich, :),1), 1:nch, 'uniformoutput',0)');
    end
end
toc

%  average across on/off gratings
RespAxTime = cellfun(@(x,y) mean(cat(3,x,y), 3), onRespAxTime, offRespAxTime, 'uniformoutput', 0);
RespAxTime_zbych = cellfun(@(x,y) mean(cat(3,x,y), 3), onRespAxTime_zbych, offRespAxTime_zbych, 'uniformoutput', 0);

fprintf('\n \t getting average responses with in specified RFtwin ... \t')
tic
% reg zscore
onResp = cellfun(@(x) nanmean(abs(x(:,et > RFtwin(1) & et < RFtwin(2))), 2), onRespAxTime, 'UniformOutput', false);
offResp = cellfun(@(x) nanmean(abs(x(:,et > RFtwin(1) & et < RFtwin(2))), 2), offRespAxTime, 'UniformOutput', false);
Resp = cellfun(@(x) nanmean(abs(x(:,et > RFtwin(1) & et < RFtwin(2))), 2), RespAxTime, 'UniformOutput', false);
% zscore bych
onResp_zbych = cellfun(@(x) nanmean(abs(x(:,et > RFtwin(1) & et < RFtwin(2))), 2), onRespAxTime_zbych, 'UniformOutput', false);
offResp_zbych = cellfun(@(x) nanmean(abs(x(:,et > RFtwin(1) & et < RFtwin(2))), 2), offRespAxTime_zbych, 'UniformOutput', false);
Resp_zbych = cellfun(@(x) nanmean(abs(x(:,et > RFtwin(1) & et < RFtwin(2))), 2), RespAxTime_zbych, 'UniformOutput', false);
toc    






medISI = median([t_info(2:end).trial_onset] - [t_info(1:end-1).trial_offset]);
if medISI > abs(diff(RFtwin))
    fprintf('\n \t getting average responses before specified RFtwin ... \t')
    tic
    % regular zscore
    onBefore= cellfun(@(x) nanmean(abs(x(:,et > -RFtwin(2) & et < RFtwin(1))), 2), onRespAxTime, 'UniformOutput', false);
    offBefore = cellfun(@(x) nanmean(abs(x(:,et > -RFtwin(2) & et < RFtwin(1))), 2), offRespAxTime, 'UniformOutput', false);
    respBefore = cellfun(@(x) nanmean(abs(x(:,et > -RFtwin(2) & et < RFtwin(1))), 2), RespAxTime, 'UniformOutput', false);
    % zscore bych
    onBefore_zbych= cellfun(@(x) nanmean(abs(x(:,et > -RFtwin(2) & et < RFtwin(1))), 2), onRespAxTime_zbych, 'UniformOutput', false);
    offBefore_zbych = cellfun(@(x) nanmean(abs(x(:,et > -RFtwin(2) & et < RFtwin(1))), 2), offRespAxTime_zbych, 'UniformOutput', false);
    respBefore_zbych = cellfun(@(x) nanmean(abs(x(:,et > -RFtwin(2) & et < RFtwin(1))), 2), RespAxTime_zbych, 'UniformOutput', false);
    toc
else
   onBefore = [];
   offBefore = [];
   respBefore  = [];
   ondRF = []; 
   offdRF = [];
   dRF = [];
   
   onBefore_zbych = [];
   offBefore_zbych = [];
   respBefore_zbych  = [];
   ondRF_zbych = []; 
   offdRF_zbych = [];
   dRF_zbych = [];
end

fprintf('\n \t converting RF tiled cells into 1xnch cells of tiled arrays... \t')
% convert cell array into 1xnch but keep tiling 
% regular zscore
onRFByCh = arrayfun(@(ich) cellfun(@(resp) resp(ich,:),onResp), 1:nch, 'uniformoutput', 0);
offRFByCh = arrayfun(@(ich) cellfun(@(resp) resp(ich,:),offResp), 1:nch, 'uniformoutput', 0);
RFByCh = arrayfun(@(ich) cellfun(@(resp) resp(ich,:),Resp), 1:nch, 'uniformoutput', 0);
% zscore bych
onRFByCh_zbych = arrayfun(@(ich) cellfun(@(resp) resp(ich,:),onResp_zbych), 1:nch, 'uniformoutput', 0);
offRFByCh_zbych = arrayfun(@(ich) cellfun(@(resp) resp(ich,:),offResp_zbych), 1:nch, 'uniformoutput', 0);
RFByCh_zbych = arrayfun(@(ich) cellfun(@(resp) resp(ich,:),Resp_zbych), 1:nch, 'uniformoutput', 0);


if ~isempty(onBefore)
    fprintf('\n \t calculating difference between baseline and stim onset... \t')
    tic
    % regular zscore
    onTemp= arrayfun(@(ich) cellfun(@(resp) resp(ich,:),onBefore), 1:nch, 'uniformoutput', 0);
    offTemp= arrayfun(@(ich) cellfun(@(resp) resp(ich,:),offBefore), 1:nch, 'uniformoutput', 0);
    rTemp = arrayfun(@(ich) cellfun(@(resp) resp(ich,:),respBefore), 1:nch, 'uniformoutput', 0);
    ondRF = cellfun(@minus, onRFByCh, onTemp, 'uniformoutput',0); 
    offdRF = cellfun(@minus, offRFByCh, offTemp, 'uniformoutput',0); 
    dRF = cellfun(@minus, RFByCh, rTemp, 'uniformoutput',0);
    
    % zscore bych
    
    onTemp_zbych= arrayfun(@(ich) cellfun(@(resp) resp(ich,:),onBefore_zbych), 1:nch, 'uniformoutput', 0);
    offTemp_zbych= arrayfun(@(ich) cellfun(@(resp) resp(ich,:),offBefore_zbych), 1:nch, 'uniformoutput', 0);
    rTemp_zbych = arrayfun(@(ich) cellfun(@(resp) resp(ich,:),respBefore_zbych), 1:nch, 'uniformoutput', 0);
    ondRF_zbych = cellfun(@minus, onRFByCh_zbych, onTemp_zbych, 'uniformoutput',0); 
    offdRF_zbych = cellfun(@minus, offRFByCh_zbych, offTemp_zbych, 'uniformoutput',0); 
    dRF_zbych = cellfun(@minus, RFByCh_zbych, rTemp_zbych, 'uniformoutput',0);
    toc
end


% get the monitor angle and elevation
[mAngle, mElev_cm] = fetchn(data.Experiments(stripKey(key(1), 'data.Experiments')), ...
    'exp_monitorangle', 'exp_monitorelevation');
% convert elevation to degrees vis angle
monitor_d = 25;
mElev = cm2deg(monitor_d, mElev_cm ,true);

% get the grating conditions (position of each grating)
[uxpos, uypos] = fetchn(data.GratingConditions(stripKey(key(1), 'data.GratingConditions')), 'grat_x_position', 'grat_y_position');
% get the regions in absolute coords (regardless of monitor angle or elevation)
sRegionX = round(unique(uxpos)) + mAngle;  % centers of tiling by squares (deg) for x
intRegionX = min(sRegionX):max(sRegionX);
sRegionY = round(unique(uypos)) + mElev;  % centers of tiling by squares (deg) for y
intRegionY = min(sRegionY):max(sRegionY);




fprintf('\n \t Interpolating RFs... \t')
tic
% regular zscore
onIntRFByCh = cellfun(@(chRF) filter2(h, chRF, 'same'), ...
    cellfun(@(chRF) interp2(sRegionX, sRegionY', chRF, intRegionX, intRegionY'),onRFByCh, ...
    'uniformoutput', 0), 'uniformoutput', 0); %#ok<NASGU>
offIntRFByCh = cellfun(@(chRF) filter2(h, chRF, 'same'), ...
    cellfun(@(chRF) interp2(sRegionX, sRegionY', chRF, intRegionX, intRegionY'),offRFByCh, ...
    'uniformoutput', 0), 'uniformoutput', 0); %#ok<NASGU>
IntRFByCh = cellfun(@(chRF) filter2(h, chRF, 'same'), ...
    cellfun(@(chRF) interp2(sRegionX, sRegionY', chRF, intRegionX, intRegionY'),RFByCh, ...
    'uniformoutput', 0), 'uniformoutput', 0); %#ok<NASGU>

% zscore bych
onIntRFByCh_zbych = cellfun(@(chRF) filter2(h, chRF, 'same'), ...
    cellfun(@(chRF) interp2(sRegionX, sRegionY', chRF, intRegionX, intRegionY'),onRFByCh_zbych, ...
    'uniformoutput', 0), 'uniformoutput', 0); %#ok<NASGU>
offIntRFByCh_zbych = cellfun(@(chRF) filter2(h, chRF, 'same'), ...
    cellfun(@(chRF) interp2(sRegionX, sRegionY', chRF, intRegionX, intRegionY'),offRFByCh_zbych, ...
    'uniformoutput', 0), 'uniformoutput', 0); %#ok<NASGU>
IntRFByCh_zbych = cellfun(@(chRF) filter2(h, chRF, 'same'), ...
    cellfun(@(chRF) interp2(sRegionX, sRegionY', chRF, intRegionX, intRegionY'),RFByCh_zbych, ...
    'uniformoutput', 0), 'uniformoutput', 0); %#ok<NASGU>





if ~isempty(onBefore)
    onIntdRF = cellfun(@(chRF) filter2(h, chRF, 'same'), ...
        cellfun(@(chRF) interp2(sRegionX, sRegionY', chRF, intRegionX, intRegionY'),ondRF , ...
        'uniformoutput', 0), 'uniformoutput', 0); %#ok<NASGU>
    offIntdRF = cellfun(@(chRF) filter2(h, chRF, 'same'), ...
        cellfun(@(chRF) interp2(sRegionX, sRegionY', chRF, intRegionX, intRegionY'),offdRF , ...
        'uniformoutput', 0), 'uniformoutput', 0); %#ok<NASGU>
    IntdRF = cellfun(@(chRF) filter2(h, chRF, 'same'), ...
        cellfun(@(chRF) interp2(sRegionX, sRegionY', chRF, intRegionX, intRegionY'),dRF , ...
        'uniformoutput', 0), 'uniformoutput', 0); %#ok<NASGU>
    
    onIntdRF_zbych = cellfun(@(chRF) filter2(h, chRF, 'same'), ...
        cellfun(@(chRF) interp2(sRegionX, sRegionY', chRF, intRegionX, intRegionY'),ondRF_zbych , ...
        'uniformoutput', 0), 'uniformoutput', 0); %#ok<NASGU>
    offIntdRF_zbych  = cellfun(@(chRF) filter2(h, chRF, 'same'), ...
        cellfun(@(chRF) interp2(sRegionX, sRegionY', chRF, intRegionX, intRegionY'),offdRF_zbych , ...
        'uniformoutput', 0), 'uniformoutput', 0); %#ok<NASGU>
    IntdRF_zbych  = cellfun(@(chRF) filter2(h, chRF, 'same'), ...
        cellfun(@(chRF) interp2(sRegionX, sRegionY', chRF, intRegionX, intRegionY'),dRF_zbych , ...
        'uniformoutput', 0), 'uniformoutput', 0); %#ok<NASGU>
else
    onIntdRF = [];
    offIntdRF = []; 
    IntdRF = []; 
    onIntdRF_zbych = [];
    offIntdRF_zbych= [];
    IntdRF_zbych= [];
end
toc

if ~isempty(matDir)
    fname = fullfile(matDir, sprintf('%s_%02d_%s_RFInfo.mat', mouse_id, key(1).series_num, expStr));
    fprintf('\n \t saving mat file %s... \t', fname)
    tic
    save(fname, ...
        'sRegionX', 'intRegionX', 'sRegionY', 'intRegionY', 'uxpos', 'uypos', ...
        'RFByCh', 'IntRFByCh', 'onRFByCh', 'onIntRFByCh', 'offRFByCh', 'offIntRFByCh', ...
        'ondRF', 'onIntdRF', 'offdRF','offIntdRF','dRF','IntdRF', ...
        'RFByCh_zbych', 'IntRFByCh_zbych', 'onRFByCh_zbych', 'onIntRFByCh_zbych', 'offRFByCh_zbych', 'offIntRFByCh_zbych', ...
        'ondRF_zbych', 'onIntdRF_zbych', 'offdRF_zbych','offIntdRF_zbych','dRF_zbych','IntdRF_zbych');
    toc
end


if plotFlag
    fprintf('\n \t plotting... godspeed... \t')
    tic
    % plot everything (all positions) combined across time 
    for rfType = {'', 'on', 'off'}
        figure('Name', sprintf('%s_%02d_%s_%sRF', mouse_id, key(1).series_num, expStr, rfType{:}));
        PlotDisplaced(et, nanmean(eval(sprintf('%sRespAxTime', rfType{:}))));
        if ~isempty(figDir)
            %print(gcf, '-depsc', fullfile(figDir, sprintf('sparseNoiseEvoked%s_%s_%02d_%s', rfType{:}, mouse_id, key(1).series_num, expStr)));
            savefig(gcf, fullfile(figDir, sprintf('sparseNoiseEvoked%s_%s_%02d_%s.fig', rfType{:}, mouse_id, key(1).series_num, expStr)));
        end
    end

    % plot the nub response for each stimulus (mean maue or min for lfp - implement plz)
    for rfType = {'', 'on','off'}
        figure('Name', sprintf('%s_%02d_%sRF%s', mouse_id, key(1).series_num, expStr, rfType{:}));
        thisRF = eval(sprintf('%sIntRFByCh', rfType{:}));
        ncols = ceil(sqrt(nch));
        nrows = ceil(nch/ncols);
        ah = nan(1, nch);
        for ich = 1 : nch
            ah(ich) = subplot(nrows,ncols,ich);
            imagesc(intRegionX, fliplr(intRegionY), thisRF {ich}); set(gca, 'YDir', 'normal');
            axis equal; colormap gray;
            set(gca, ...
                'XTick', [sRegionX(1) sRegionX(end)], ...
                'YTick', [sRegionY(1) sRegionY(end)]);
        end
        allAxLims = cell2mat(arrayfun(@caxis, ah, 'UniformOutput', false)');
        minLim = min(allAxLims(:,1));
        maxLim = max(allAxLims(:,2));
        arrayfun(@(x)(caxis(x, [minLim, maxLim])), ah, 'UniformOutput', false)
        if ~isempty(figDir)
            %print(gcf, '-depsc', fullfile(figDir, sprintf('sparseNoiseMap%s_%s_%02d_%s', rfType{:}, mouse_id, key(1).series_num, expStr)));
            savefig(gcf, fullfile(figDir, sprintf('sparseNoiseMap%s_%s_%02d_%s.fig', rfType{:}, mouse_id, key(1).series_num, expStr)));
        end
    end
    toc
end










