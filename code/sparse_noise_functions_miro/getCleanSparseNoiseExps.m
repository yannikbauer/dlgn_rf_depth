function cleanExps = getCleanSparseNoiseExps(varargin)
% cleanExps = getCleanSparseNoiseExps(varargin)
% cleanExps = getCleanSparseNoiseExps returns everything
% cleanExps = getCleanSparseNoiseExps(series_key) only returns the series
% cleanExps = getCleanSparseNoiseExps({mouse_id, series_num}) same as previous
% cleanExps = getCleanSparseNoiseExps(mouse_key) only returns the mouse
% cleanExps = getCleanSparseNoiseExps('lgn') only returns lgn series
% cleanExps = getCleanSparseNoiseExps(anyOtherString) only returns v1 series
% cleanExps = getCleanSparseNoiseExps

narginchk(0,1)

%%%%%%%%%%%%%%%%%
%%%% V1,ChR2 %%%%
%%%%%%%%%%%%%%%%%

% Ntsr1-Cre_0080
allCleanExps = struct('mouse_id', 'Ntsr1-Cre_0080', ...
    'series_num', 2, 'cleanExps', 6 ,'badExps',[1,4], 'isLgn', false, 'badChans', [16:19,29:32]);

% Ntsr1-Cre_0082
allCleanExps(end+1) = struct('mouse_id', 'Ntsr1-Cre_0082', ...
    'series_num', 2, 'cleanExps', [4,14] ,'badExps',2,'isLgn', false, 'badChans', [1:2, 16:22, 26:32]);

%%%%%%%%%%%%%%%%%%
%%%% lgn,ChR2 %%%%
%%%%%%%%%%%%%%%%%%











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% add mouse_counter in the end  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mcs = cellfun(@(m) ...
    fetch1(data.Mice(sprintf('mouse_id like "%%%s%%"', m)), 'mouse_counter'),...
    {allCleanExps.mouse_id}, 'uniformoutput',0); 

[allCleanExps.mouse_counter] = deal(mcs{:});
clear mcs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% deterimining output %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    myLog = true(size(allCleanExps));
elseif isstruct(varargin{1})
    searchKey = varargin{1};
    if isfield(searchKey, 'series_num')
        myLog = arrayfun(@(s) isequal(s, searchKey), stripKey(allCleanExps, 'data.Series'));
    else
        myLog  = arrayfun(@(s) isequal(s, searchKey), stripKey(allCleanExps, 'data.Mice'));
    end
elseif iscell(varargin{1})
    searchKey = varargin{1};
    myLog = strcmpi({allCleanExps.mouse_id}, searchKey{1}) &  [allCleanExps.series_num]==searchKey{2};
elseif ischar(varargin{1}) 
    myLog  = [allCleanExps.isLgn] == strcmpi(varargin{1},'lgn');
else
    error('undefined input')
end

cleanExps = allCleanExps(myLog);