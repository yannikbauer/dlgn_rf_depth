% script to keep track of cleaning sparseNoiseExperiments
% will use call this to make allCleanExps

%%%%%%%%%%%%%%%%%
%%%% V1,ChR2 %%%%
%%%%%%%%%%%%%%%%%
allCleanExps = struct('mouse_id', 'Ntsr1-Cre_0080', ...
    'series_num', 2, 'cleanExps', 6 ,'badExps',[1,4], 'lgn', false);

allCleanExps(end+1) = struct('mouse_id', 'Ntsr1-Cre_0082', ...
    'series_num', 2, 'cleanExps', [4,14] ,'badExps',2,'lgn', false);



% add mouse counter in the end
mcs = cellfun(@(m) ...
    fetch1(data.Mice(sprintf('mouse_id like "%%%s%%"', m)), 'mouse_counter'),...
    {allCleanExps.mouse_id}, 'uniformoutput',0); 

[allCleanExps.mouse_counter] = deal(mcs{:});
clear mcs