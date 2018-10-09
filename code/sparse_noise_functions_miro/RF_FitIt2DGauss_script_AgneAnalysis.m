%%% ignore (not true) 
% diff from original: 
% WHEN DEFINING "myGaussian", 
% scaletoOne OPTION FOR TwoDGaussian.m 
% IS TRUE 
% may make more sense in summing (double check the sum code) 
%%%%%

% Fits with Fitit (how we usually do it)
close all; 
clear all; startup;
addpath(genpath('~/Dropbox/clean/locomotion_layer/animals/'));
%rmpath(genpath('~/Dropbox/clean/locomotion_layer/animals/agneReAnalysis/sparseNoiseFunctions/obsolete/2dgaussian301'));

onlyPlotFits = false;
%muaDir = '/Users/yoknapatawpha/Desktop/SparseNoise_LGN/sparseNoiseMats_MUAE';
% muaDir = '/Volumes/busse_lab/users/sinem/SparseNoise_LGN/sparseNoiseMats_MUAE';
% muaDir = '/Volumes/busse_lab/users/sinem/SparseNoise_V1/sparseNoiseMats_MUAE';
muaDir = '/Volumes/busse_lab/users/sinem/SparseNoise_agneSS/sparseNoiseMats_MUAE';
% options for selecting good fits while plotting
rsq_co = .4; %0.42 ; % for plotting 
minRadius = 1;  % for sigmax, sigmay
maxRadius = 16; % for sigmax, sigmay

if ~onlyPlotFits
    plotNansOnly = false;
else
    plotNansOnly = false;
    % for visualization
    setFigShit
    zPos = [33, 16, 0];
    rfPos = [0 30 60];
    subplotsAreTight = true;
    scrnfigw = 24;
    scrnfigh = 16;
    
    %%
    %close all;
%     set(groot, 'Units', figPars.Units)
%     scrnsz = get( groot, 'Screensize');
%     scrnps = get( groot, 'MonitorPosition');
%     
%     scrnfigw =20;
%     scrnfigh = 16;
%     nzpos = 2;
%     nrfpos = 3;
%     zPos = linspace(scrnsz(4)+2-scrnfigh,0,nzpos)+scrnps(2) 
%     rfPos = linspace(0,max(scrnsz(3)-scrnfigw, scrnfigw*(nrfpos -1)), nrfpos )
%     subplotsAreTight = true;
%     zfigidx = 1;
    
%     arrayfun(@(i) figure(figPars, 'Name', num2str(i),'position', [rfPos(i) zPos(i), scrnfigw,scrnfigh], ...
%         'paperposition', [rfPos(i) zPos(i), figsz,scrnfigh]), 1:nfigs)
    %%
end

% TILT ALLOWED  (ROTATION)
noTilt = false;
% instead of rotating, use chol factor to build cov matrix
cholFactor = false;
% make dirs
if noTilt
    fitDir = fullfile(muaDir, 'fitGauss2D');
elseif cholFactor
    fitDir = fullfile(muaDir, 'fitGauss2D_chol');
    %zfigidx = 1;
    zfigidx = [];
else
    fitDir = fullfile(muaDir, 'fitGauss2D_tilted');
    zfigidx = [];
    
end

if ~isdir(fitDir)
    mkdir(fitDir)
end


% rfTypes
rfType = {'', 'on', 'off'};

%
z_cos = 20; [8 10 20]; %[10 8 20]; %[8,10,20]; *8,10 already running (mcmc, positive)
changeZpos = false; % for plotting on the monitor 
for z_co = z_cos 

    if isempty(zfigidx) || changeZpos
        changeZpos = true;
        zfigidx = find(z_cos == z_co);
    end
    
    zFitDir = fullfile(fitDir, sprintf('z_co%d',z_co));
    if ~isdir(zFitDir)
        mkdir(zFitDir);
    end
    
    
    zDir = fullfile(muaDir, sprintf('z_co%d',z_co));
    tmp = dir(fullfile(zDir, '*_RFInfo*'));
    [~, tmp1] = regexp({tmp.name}, '([\S]+)\_RFInfo.mat', 'match', 'tokens', 'once');
    
    basenames = [tmp1{:}];
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%% for rerunning Ntrs1-Cre_0080s07 %%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     baseLog = cellfun(@isempty, strfind(basenames, 'Ntsr1-Cre_0080_07'));
%     basenames = basenames(~baseLog);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % inline function for calling TwoDGaussian
    if exist('myGaussian', 'var'); clear myGaussian; end;
    if cholFactor
        myGaussian = @(pars, coords) TwoDGaussian(pars, coords, true, true);
    else
        myGaussian = @(pars, coords) TwoDGaussian(pars, coords, true, false);
    end
    
    
    
    for ibase = 1:numel(basenames)
        
        RF = load(fullfile(zDir, [basenames{ibase},'_RFInfo.mat']));
        nch = numel(RF.IntRFByCh);
        
        if ~onlyPlotFits
            
            rng(1) % for reproducibility of fits
            
            
            %if size(RF.intRegionY,1) ~= 1; RF.intRegionY = RF.intRegionY'; end
            [xi,yi] = meshgrid(RF.intRegionX,RF.intRegionY);
            myFit = struct('IntRFByCh',  cell(1,nch), 'onIntRFByCh', cell(1,nch),'offIntRFByCh', cell(1,nch), 'ch', num2cell(1:nch));
            
            for ich = 1:nch
                for irf = 1:numel(rfType)
                    itype = 'IntRFByCh';
                    if ~strcmp(rfType{irf}, 'rf')
                        itype = [rfType{irf} itype]; %#ok<AGROW>
                    end
                    
                    
                    % original image
                    z = flipud(RF.(itype){ich});
                    
                    
                    % initialize structure with fit info
                    fitStruct = struct('pars', {[]},'fitImg', {[]}, ...
                        'r2',{nan}, 'sse', {nan}, 'sst', {sum((z(:)-mean(z(:))).^2)},...
                        'x_fit', {[]}, 'y_fit', {[]}, 'x0', {nan}, 'y0',{nan}, 'C', {[]});

                    
                    if any(isnan(z(:)))
                        figure('Name', [basenames{ibase} sprintf('z_co%d',z_co)])
                        imagesc(RF.intRegionX,RF.intRegionY,z);
                        myFit(1,ich).(itype) = fitStruct;
                        continue
                    elseif plotNansOnly
                        continue
                    end
                    
                    
                    % starting parameters
                    % center: find max
                    [max_yi,max_xi] = find(z==max(z(:)));
                    x0 = RF.intRegionX(max_xi);
                    y0 = RF.intRegionY(max_yi);
                    % offset
                    offset_0 = prctile(z(:), 10);
                    % amplitude
                    amp_0 = max(z(:))-offset_0;
                    % tilt
                     if noTilt 
                        maxAngle = 0;
                        theta_0 = 0;
                     elseif cholFactor
                         % angle doesnt mean anything;
                         theta_0 = 0;
                         maxAngle = range(RF.intRegionX);
                     else
                        maxAngle = 180;
                        theta_0= 90;
                     end

                    % spread
                    sigmax_0 = 5;
                    sigmay_0 = 5;
                    % starting pars as vector
                    pars_0 = [x0, y0, sigmax_0, sigmay_0, amp_0,offset_0,theta_0];
                    
                    % lower bound for parameters
                    center_out = 10; % how much the center can be out of sparse noise region
                    pars_lb = [min(RF.intRegionX)-center_out, min(RF.intRegionY)-center_out, ...
                                0.5, 0.5, 0, 0, 0];
                    
                    % upper bound for parameters
                    pars_ub = [max(RF.intRegionX)+center_out, max(RF.intRegionY)+center_out ,...
                                range(RF.intRegionX), range(RF.intRegionY), max(z(:)),2*max(z(:)), maxAngle];
                    % options 
                    %[display parameter, tol x, tol dist, #start pts, niter, lastp is weights]
                    fitit_opts = [1 1e-4 1e-4 10 1e4 0];
                    [~, pars] = fitit(myGaussian, z, ...
                                        pars_lb, pars_0, pars_ub, fitit_opts,...
                                        [{RF.intRegionX} {RF.intRegionY}]);
                    [fit,C] = myGaussian(pars,[{RF.intRegionX}, {RF.intRegionY}]);
                    SSE = sum((z(:)-fit(:)).^2);
                    r2 = 1-(SSE/fitStruct.sst);
                    
                    % get elipse
                    % from: http://www.cs.ubc.ca/~murphyk/Teaching/CS340-Fall06/reading/gaussians.pdf
                    mu = pars(1:2)';
                    [x_elipse, y_elipse] = getGaussEllipse(mu,C, 0.95);

                    
                    fitStruct.pars = pars;
                    fitStruct.fitImg = fit;
                    fitStruct.r2 = r2;
                    fitStruct.sse = SSE;
                    fitStruct.x_fit = x_elipse;
                    fitStruct.y_fit = y_elipse;
                    
                    fitStruct.x0 = pars(1);
                    fitStruct.y0 = pars(2);
                    fitStruct.C = C;
                    
                    % plot for testing
%                     figure; h(1)= subplot(1,2,1); imagesc(RF.intRegionX, RF.intRegionY, z); h(2) = subplot(1,2,2); imagesc(RF.intRegionX, RF.intRegionY,fit);
%                     set(h, 'clim', [min(cat(1,z(:), fit(:))) max(cat(1,z(:), fit(:)))], 'ydir', 'normal')
%                     hold(h(2), 'on')
%                     plot(h(2), x_elp, y_elp,'r')
                    % 
                    myFit(1,ich).(itype) = fitStruct;
                end
            end
            fname = fullfile(zFitDir, [basenames{ibase},'_fits.mat']);
            save(fname, 'myFit')
            
        elseif onlyPlotFits  % only plot autoFits (based off directory)
            fname = fullfile(zFitDir, [basenames{ibase},'_fits.mat']);
            if ~exist(fname, 'file')
                warning('cant find %s', fname);
                continue
            end
            tmp = load(fname);
            myFit = tmp.myFit; % gaussian fits
            %%
            %get the RFs
            [xi, yi] = meshgrid(RF.intRegionX, RF.intRegionY);                    
            for irf = 1:numel(rfType)
                itype = rfType{irf};
                % each type organized by ch
                thisRF = [itype 'IntRFByCh'];
                
                % stupid bug while making mats ... fixed code now but mats need
                % to be rerun
                badLog = ~arrayfun(@(ch) isfield(ch.(thisRF), 'C'), myFit);
                if any(badLog) % nan rfs dont have a C
                    fitStruct = struct('pars', {[]},'fitImg', {[]}, ...
                        'r2',{nan}, 'sse', {nan}, 'sst', {nan},...
                        'x_fit', {[]}, 'y_fit', {[]}, 'x0', {nan}, 'y0',{nan}, 'C', {[]});
                    
                    [myFit(badLog).(thisRF)] = deal(fitStruct);
                end
                
                
                fitRF = [myFit.(thisRF)]; % gauss
                if numel(fitRF) ~= nch 
                      error('deal with me. fitRF has missing channels')
                end
                
                figure(figPars, 'Name', sprintf('z_co=%d %s %s',z_co, basenames{ibase}, thisRF), ...
                    'position', [rfPos(irf) zPos(zfigidx), scrnfigw,scrnfigh])
                h = plotRFbyChStruct(RF.(thisRF), RF.sRegionX, RF.sRegionY,RF.intRegionX,RF.intRegionY,subplotsAreTight);
                for ich = 1:nch
                    if isnan(fitRF(ich).r2)
                        continue % if r2 is nan
                    else
                        % axes: hold on to draw elipse
                        hold(h(ich), 'on')
                    end
                    
                    
                    % determine if good fit
                    R2 = fitRF(ich).r2;
                    
%                     goodFit = R2 >rsq_co &&  prod(fitRF(ich).pars(3:4)) < maxRadius^2 ...
%                                         && prod(fitRF(ich).pars(3:4)) > minRadius^2;
                                    
                    goodFit = R2 >rsq_co;

                    if goodFit
                        plot(h(ich),fitRF(ich).x_fit, fitRF(ich).y_fit,'r-')
                        plot(h(ich), fitRF(ich).x0,fitRF(ich).y0 , 'r+');
                    elseif rsq_co > 0.5 && R2 > 0.5
                        plot(h(ich),fitRF(ich).x_fit, fitRF(ich).y_fit,'g-')
                        plot(h(ich), fitRF(ich).x0,fitRF(ich).y0 , 'g+');
                    elseif rsq_co > 0.4 && R2 > 0.4
                        plot(h(ich),fitRF(ich).x_fit, fitRF(ich).y_fit,'c-')
                        plot(h(ich), fitRF(ich).x0,fitRF(ich).y0 , 'c+');
                    else
                        plot(h(ich),fitRF(ich).x_fit, fitRF(ich).y_fit,'b-')
                        plot(h(ich), fitRF(ich).x0,fitRF(ich).y0 , 'b+');
                    end
 
                    if subplotsAreTight
                        mytitle = cellfun(@(s,val) sprintf(s,val), {'%d', '%d'},...
                            {ich, round(100*R2)}, 'uniformoutput',0);
                    else
                        mytitle = sprintf('%d, R2=%0.2f',  ich, R2);
                    end
                    
                    if goodFit                        
                        title(h(ich),mytitle, 'color', 'r');
                    elseif R2 > rsq_co
                        title(h(ich),mytitle, 'color', 'b');
                    else
                        title(h(ich),mytitle);
                    end
                end
            end
             
        end %
        %pause; close all;
    end % ibase
    
end % z_co

