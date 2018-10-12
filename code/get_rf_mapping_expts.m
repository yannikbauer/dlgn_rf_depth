%% dLGN RF depth analysis
% Gets RF mapping experiments (using sparse noise stimulus) for units used
% in Roman_Roson_el_al and writes track file used for later export of
% filtered MUA recordings.

%% Project path
cd(fileparts(which(mfilename))); % changes working dir to current file
addpath(genpath('../../')); % add entire project path, incl. general code used across analyses

%% Check DataJoint connection
if ~exist('DJ_DIRS')
    startup_lmu;
end

%% Parameters
saveDir = '../data/';

%% GET KEYS USED FOR UNITS IN OTHER ANALYSES OF ROMAN_ROSON_ET_AL_2018
corr_p      = 0.001;         % 0.001; corr_p represents the wilcoxon ransum correlation for between/withing segments
qi          = 0.05;          % 0.05; Philipps quality index
cluster_qi  = 3;             % 3; cluster quality index of likely single- vs multiunit
r           = 0.65;          % 0.7; regression index
type        = 3;

unitList = getAllUnits(corr_p, qi, cluster_qi, r, type);
% save(fullfile(saveDir, 'dlgn_unit_list.mat'), 'unitList');

%% Get corresponding sparse noise experiments

% Get experiment numbers of sparse noise experiments
sn_exp = fetch(data.Experiments('exp_name LIKE "AsparseNoise%"') & stripKey(unitList, 'data.Series'));

% Add corresponding mouse IDs/name and date
sn_exp = fetch(data.Mice * data.Experiments(sn_exp) * data.Series,  'mouse_id', 'series_date');

% Sort by date
[~, idx] = sort({sn_exp.series_date});
sn_exp = sn_exp(idx);

% Check that all data use same electrode
array_num = fetch(data.Series(sn_exp), 'array_num');
array_model = fetchn(data.Arrays(array_num), 'array_model');
if length(unique(array_model)) ~= 1
    error("Experiments use different array models!") 
end


%% Get file names of sparse noise .ns6 files
% e.g. Ntsr1-Cre_0027_s06_20140408_01.ns6

for i = 1:length(sn_exp)
    % [animalFolder, seriesFolder, unsortedDataFolder, sortedDataFolder] = getEphysDataFolder(animalName, seriesNum, unsortedDataSuffix);
    [~, ~, unsortedDir, dir] = getEphysDataFolder(sn_exp(i).mouse_id, sn_exp(i).series_num); % get data folder
    dir_parts = strsplit(dir, '/'); % extract last dir part for file id
    fid = strcat(dir_parts{end}, '_', sprintf('%02d', sn_exp(i).exp_num), '.ns6'); % construct file id
    
    sn_exp(i).fid = fid; % insert into struct
    
    % Check if file exists
    if ~exist(fullfile(unsortedDir,fid), 'file')
        error("File does not exist.")
    end
end

save(fullfile(saveDir, 'sn_exp.mat'), 'sn_exp');

%% Write file names to text file (.track)

fileID = fopen(fullfile(saveDir, 'sn_file_list.track'), 'w');
formatSpec = '%s\n';

for i = 1:length(sn_exp)
    fprintf(fileID, formatSpec, sn_exp(i).fid);
end

fclose(fileID);

