function [unitList] = getAllUnits(corr_p,qi,cluster_qi,r,type)
% unitList = getAllUnits(corr_p,qi,cluster_qi,type) The script creates a
% table of all units that satisfy specific quality criteria in their 
% response to the chirp stimulus.
%
%   corr_p      = corr_p represents the wilcoxon ransum correlation for between/withing segments
%   qi          = Philipps quality index
%   cluster_qi  = cluster quality index of likely single - vs multiunit
%   r           = regression index
%   type        = determines the uni types: 
%
%                 1 = chirp units
%                 2 = SFTFORI units
%                 3 = both  
%
%   [unitList] = getAllUnits(...) returs a list of units
%
if nargin < 4
    corr_p = 0.001;
    qi = 0.05;
    cluster_qi = 3;
    r = 0.65;
    type = 3;    
end

fprintf('\nFUNCTION: GET UNITS\n')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%  CHIRP KEYS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if type == 1 || type == 3
    fprintf('... Getting chirp units\n');
    
    % get the chirp keys but exclude all experiments before '2014-03-14' and all ntsr-mice 
    chirp_keys = (data.Series('series_date > "2014-03-14"') & ...
        data.Experiments('exp_name LIKE "%chirp%"')); %- data.Mice('mouse_id LIKE "%ntsr%"');
    
    % Exclude all experiments before '2014-03-14'
    date_keys = fetch(data.Series(chirp_keys),'series_date');
    mask_date = datenum(fetchn(data.Series(chirp_keys),'series_date')) > datenum('2014-03-14', 'yyyy-mm-dd');
    
    % Contains all valid chirp keys
    chirp_keys = date_keys(mask_date);
    chirp_keys = rmfield(chirp_keys,'series_date');
    
    % Create a table of units that respond satisfactorily to chirp stimulus
    chirpUnits = fetch(miro.ChirpQuality(chirp_keys) ...
        & data.ClusterInfo(sprintf('quality <= %d',cluster_qi)) ...        
        & miro.ChirpQuality(sprintf('corr_p <= %f',corr_p)) ...
        & miro.ChirpQuality(sprintf('berens_qi >= %f',qi)), ...     
        'corr_p', 'berens_qi');
    
    % return chirp units for type 1
    if type == 1
        unitList = chirpUnits;
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%  SFTFOri KEYS  %% %%%%%%%%%%%%%%%%%%%%%%%%%%

if type == 2 || type == 3
    
    % COMPUTE VISULA RESPONSIVENESS %
    fprintf('... Getting SFTFOri units\n');
    
    % x Number of points have to be over the threshold
    nVis = 10; % 10% -> 10
    allUnits = fetch(food.SFTFOriUnitStability(sprintf('r>=%f',r)) & data.ClusterInfo(sprintf('quality <= %f',cluster_qi)), 'slope');
    unit_mask = zeros(numel(allUnits), 1);
    
    for iunit = 1:numel(allUnits)
        
        stimInfo = fetch(data.StimInfo(allUnits(iunit)) & fetch(data.SFTFOriTuning(allUnits(iunit))), '*'); % contains both SFTFORIexps
        
        % check for 2 experiments
        if numel(stimInfo) == 1
            continue
        end
        
        % Go through each experiment separately
        visResExp = zeros(2,1);
        for iexp = 1:2
            
            key = allUnits(iunit);
            key.exp_num = stimInfo(iexp).exp_num;
            
            c = fetch(data.ConditionSpikes(key), 'grat_num', 'cond_rate', 'cond_sem');
            activeLog = ismember([c.grat_num], [stimInfo.num_active_grats]);
            blankLog  = ismember([c.grat_num], [stimInfo.num_blank_grats]);
            
            % positive response || negative response
            if nnz([c(activeLog).cond_rate] - 2.58*[c(activeLog).cond_sem] > mean([c(blankLog).cond_rate])) > nVis || ...
                    nnz([c(activeLog).cond_rate] + 2.58*[c(activeLog).cond_sem] < mean([c(blankLog).cond_rate])) > nVis
                visResExp(iexp) = true;
            end
        end
        
        if(visResExp(1) == 1 && visResExp(2) == 1)
            unit_mask(iunit) = true;
        end
    end
    
    allUnits(unit_mask == 0) = [];
    
    
    
    % FILTER FIRING RATES
    threshold = 2;
    allUnits(([allUnits.slope] < 1/threshold) & ([allUnits.slope] > 1*threshold)) = [];
    
    % return SFTFOri units for type 2
    if type == 2
        unitList = allUnits;
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%  JOIN BOTH UNIT LISTS  %%%%%%%%%%%%%%%%%%%%%%%%%%

if type == 3
    
    fprintf('... Creating joint unit list\n')
    
    % bring structures to the same baseline
    % adds chirp exp_num
    chirpUnits = stripKey(chirpUnits, 'data.Units');
    allUnits = stripKey(allUnits, 'data.Units');
    allUnits = fetch(data.Units(allUnits) & fetch(data.Experiments(allUnits) & 'exp_name LIKE "%Chirp%"'));
    
    % same units from chirpUntis will be removed in allUnits
    allUnits = fetch(data.Units(allUnits) - data.Units(chirpUnits));
    
    % add both structures together
    unitList = [allUnits' chirpUnits']';
    
end

unitList = fetch(data.Units(unitList));
end
























