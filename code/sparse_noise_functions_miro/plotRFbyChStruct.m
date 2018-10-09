function ah =  plotRFbyChStruct(thisRF, sRegionX, sRegionY,intRegionX,intRegionY,tightFlag, RFcenters)
% plotRFbyChStruct(thisRF, sRegionX, sRegionY,intRegionX,intRegionY,tightFlag, RFcenters, Fitxy)

if nargin < 7 || isempty(RFcenters)
    plotCenters = false;
else
    plotCenters = true;
end

if nargin < 6 || isempty(tightFlag)
    tightFlag = false;
end

nch = numel(thisRF);
ncols = ceil(sqrt(nch));
nrows = ceil(nch/ncols);
ah = nan(1, nch);
for ich = 1 : nch
    
    if tightFlag
        ah(ich) = subplottight(nrows,ncols,ich);
    else
        ah(ich) = subplot(nrows,ncols,ich);
    end
    
    imagesc(intRegionX, fliplr(intRegionY), thisRF {ich}); set(gca, 'YDir', 'normal');
    axis equal; colormap gray;
    th = title(num2str(ich));
    % Change title position
    if tightFlag
        v= axis;
        set(th, 'HorizontalAlignment', 'left', 'VerticalAlignment','top', 'Position', [v(1) v(4)]);
    end
    set(gca, ...
        'XTick', [ceil(sRegionX(1)) floor(sRegionX(end))], ...% 'XTickLabel',[], ...
        'YTick', [ceil(sRegionY(1)) floor(sRegionY(end))]);  % 'YTickLabel',[]);
    
    if plotCenters && all(~isnan(RFcenters(ich,:)))
        plot(ah(ich), RFcenters(ich,1), RFcenters(ich,2), 'marker', 'r+', 'markersize', 5);
    end
end
allAxLims = cell2mat(arrayfun(@caxis, ah, 'UniformOutput', false)');
minLim = min(allAxLims(:,1));
maxLim = max(allAxLims(:,2));
arrayfun(@(x)(caxis(x, [minLim, maxLim])), ah, 'UniformOutput', false)