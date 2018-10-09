% Fits with Fitit (how we usually do it)
%close all; clear all; startup;
addpath(genpath('~/Dropbox/clean/locomotion_layer/animals/'));
rmpath(genpath('~/Dropbox/clean/locomotion_layer/animals/agneReAnalysis/sparseNoiseFunctions/2dgaussian301'));

muaDir = '/Users/yoknapatawpha/Desktop/SparseNoise_LGN/sparseNoiseMats_MUAE';
%muaDir = '/Volumes/busse_lab/users/sinem/SparseNoise_LGN/sparseNoiseMats_MUAE';

onlyPlotFits = true; % doesnt compute fits if set to true
% save RFcenters (for onlyPlotFits = true)

rsqCo = 0.42 ; % for plotting

% options for computing RFs: errorbar
if ~onlyPlotFits
    plotNansOnly = false;
else
    setFigShit
    zPos = [33, 16, 0];
    rfPos = [0 30 60];
end

% make dirs
fitDir = fullfile(muaDir, 'fitEllipse');
if ~isdir(fitDir)
    mkdir(fitDir)
end


% rfTypes
rfType = {'', 'on', 'off'};

%

for z_co =10 % [10 8 20]; %[8,10,20]; *8,10 already running (mcmc, positive)
    
    zFitDir = fullfile(fitDir, sprintf('z_co%d',z_co));
    if ~isdir(zFitDir)
        mkdir(zFitDir);
    end
    
    
    zDir = fullfile(muaDir, sprintf('z_co%d',z_co));
    tmp = dir(fullfile(zDir, '*_RFInfo*'));
    [~, tmp1] = regexp({tmp.name}, '([\S]+)\_RFInfo.mat', 'match', 'tokens', 'once');
    basenames = [tmp1{:}];
    
    for ibase = 1:numel(basenames)
        RF = load(fullfile(zDir, [basenames{ibase},'_RFInfo.mat']));
        nch = numel(RF.IntRFByCh);
        
        if ~onlyPlotFits
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
                        'x_ellipse', {[]}, 'y_ellipse', {[]}, 'x0', {nan}, 'y0',{nan});
                    
                    
                    
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
                    theta_0 = pi/2;
                    % spread
                    sigmax_0 = 5;
                    sigmay_0 = 5;
                    % starting pars as vector
                    pars_0 = [x0, y0, theta_0, amp_0, sigmax_0, sigmay_0, offset_0];
                    
                    % lower bound for parameters
                    center_out = 10; % how much the center can be out of sparse noise region
                    pars_lb = [min(RF.intRegionX)-center_out, min(RF.intRegionY)-center_out, ...
                                0, 0, 0, 0, 0];
                    % upper bound for parameters
                    pars_ub = [max(RF.intRegionX)+center_out, max(RF.intRegionY)+center_out ,...
                                pi, 2*max(z(:)),  range(RF.intRegionX), range(RF.intRegionY), max(z(:))];
                    % options 
                    %[display parameter, tol x, tol dist, #start pts, niter, lastp is weights]
                    fitit_opts = [1 1e-4 1e-4 10 1e4 0];
                    [~, pars] = fitit('twoDEllipseWithOffset', z, ...
                                        pars_lb, pars_0, pars_ub, fitit_opts,...
                                        [{RF.intRegionX} {RF.intRegionY}]);
                    fit = twoDEllipseWithOffset(pars,[{RF.intRegionX} {RF.intRegionY}]);
                    SSE = sum((z(:)-fit(:)).^2);
                    r2 = 1-(SSE/fitStruct.sst);
                    
                    % elipse needs angle in degrees
                    angleDeg = pars(3)*180/pi;
                    angleDeg  = degdiff(180, angleDeg, 180);
                    [x_elp, y_elp] = calculateEllipse(pars(1), pars(2), pars(5), pars(6), angleDeg);
                    
                    fitStruct.pars = pars;
                    fitStruct.fitImg = fit;
                    fitStruct.r2 = r2;
                    fitStruct.sse = SSE;
                    fitStruct.x_ellipse = x_elp;
                    fitStruct.y_ellipse = y_elp;
                    fitStruct.x0 = pars(1);
                    fitStruct.y0 = pars(2);
                    
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
            myFit = tmp.myFit; % eliptical fits
            %%
            %get the RFs
            [xi, yi] = meshgrid(RF.intRegionX, RF.intRegionY);                    
            for irf = 1:numel(rfType)
                itype = rfType{irf};
                % each type organized by ch
                thisRF = [itype 'IntRFByCh'];
                fitRF = [myFit.(thisRF)]; % gauss
                if numel(fitRF) ~= nch 
                      error('deal with me')
                end
                figure(figPars, 'Name', sprintf('z_co=%d %s %s',z_co, basenames{ibase}, thisRF), ...
                    'position', [rfPos(irf) zPos(2), 30,30])
                h = plotRFbyChStruct(RF.(thisRF), RF.sRegionX, RF.sRegionY,RF.intRegionX,RF.intRegionY);
                for ich = 1:nch
                    if isnan(fitRF(ich).r2)
                        continue % if r2 is nan
                    else
                        % axes: hold on to draw elipse
                        hold(h(ich), 'on')
                    end

                    % determine if good fit
                    R2 = fitRF(ich).r2;
                    goodFit = R2 >rsqCo && prod(fitRF(ich).pars(5:6)) < 16^2  && prod(fitRF(ich).pars(5:6)) > 1;

                    if goodFit
                        plot(h(ich),fitRF(ich).x_ellipse, fitRF(ich).y_ellipse,'r-')
                        plot(h(ich), fitRF(ich).x0,fitRF(ich).y0 , 'r+');
                    else
                        plot(h(ich),fitRF(ich).x_ellipse, fitRF(ich).y_ellipse,'b-')
                        plot(h(ich), fitRF(ich).x0,fitRF(ich).y0 , 'b+');
                    end
 
                    mytitle = sprintf('%d, R2=%0.2f',  ich, (R2));
                    if goodFit                        
                        title(h(ich),mytitle, 'color', 'r');
                    elseif R2 > rsqCo
                        title(h(ich),mytitle, 'color', 'b');
                    else
                        title(h(ich),mytitle);
                    end
                end
            end
            pause
        end %
    end % ibase
end % z_co

