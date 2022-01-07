clear all; close all; clc
addpath('functions')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% title : mainSproutQuant.m
% Description : analyses the neovascularization sprouts density and length 
%               as a function of their localization.
% Author : Nicolas Desjardins-Lecavalier
% Edition date : 26 novembre 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


prompt = {'Enter the image folder','Enter the image resolution (microns/pixel)','Enter the gaussian filter size (pixels)','Enter the estimated sprout size (pixels)'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'C:\Users\Nicolas\Dropbox (Biophotonics)\In Vivo CLaP\Sample images for sprout quant','','15','5'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

imFolder = answer{1};
imRes = str2num(answer{2});
gaussFiltSize = str2num(answer{3});
sproutSize = str2num(answer{4});

resultFolder = fullfile(imFolder,'results');

if not(exist(resultFolder))
    mkdir(resultFolder)
end

imFiles = dir(fullfile(imFolder,'*.tif'));
resultTable = [];

% for idx = 1:numel(imFiles)
for idx =1:2
    t = imFiles(idx).name;
    V = tiffreadVolume(fullfile(imFolder,t));
    I = mat2gray(max(V,[],3)); % max projection
    Ihisteq = adapthisteq(I); %
    Ihighpass = imgaussfilt(Ihisteq,1)-imboxfilt(Ihisteq,sproutSize*2+1) ;
    
    % finds skeleton
    mask = imbinarize(Ihighpass,'adapt');
    mask = imerode(mask,strel('disk',1));
    mask = bwareaopen(mask,10);
    skeleton = bwskel(mask,'minbranchlength',11);
    [skeletonLbl,nSkelLbl] = bwlabel(skeleton);
    
    % finds contours
    roughMask = imbinarize(imgaussfilt(Ihisteq,gaussFiltSize));
    roughMask = imerode(roughMask,strel('disk',round(gaussFiltSize/2)));
    roughMask = imfill(roughMask,'holes');
    roughContour = bwperim(roughMask);
    contourDist1 = imdilate(roughMask,strel('disk',10));
    contourDist1 = bwperim(contourDist1);
    [contourDist1Lbl, nContLbl] = bwlabel(contourDist1);
    
    % finds intersections (x0,y0)
    [X,Y] = meshgrid(1:size(I,1),1:size(I,2));
    
    xcontour = [];
    ycontour = [];
    for iLbl = 1: nContLbl
        xLbl = X(contourDist1Lbl==iLbl);
        yLbl = Y(contourDist1Lbl==iLbl);
        
        [xLblSorted,yLblSorted] = untangle(xLbl,yLbl);
        
        xcontour = [xcontour;xLblSorted];
        xcontour = [xcontour;NaN];
        ycontour = [ycontour;yLblSorted];
        ycontour = [ycontour;NaN];
    end
    
    xskeleton = [];
    yskeleton = [];
    for iLbl = 1: nSkelLbl
        xLbl = X(skeletonLbl==iLbl);
        yLbl = Y(skeletonLbl==iLbl);
        
        [xLblSorted,yLblSorted] = untangle(xLbl,yLbl);
        
        xskeleton = [xskeleton;xLblSorted];
        xskeleton = [xskeleton;NaN];
        yskeleton = [yskeleton;yLblSorted];
        yskeleton = [yskeleton;NaN];
    end
    
    [x0,y0,iout,jout] = intersections(xcontour,ycontour,xskeleton,yskeleton,0);
    groups = dbscan([x0,y0],10,1);
    
    % identifies what is a sprout
    sproutstmp = skeleton & not(roughMask);
    sprouts = false(size(sproutstmp));
    [sproutsLbltmp,nsproutsLbl] = bwlabel(sproutstmp);
    
    sproutLocs = zeros(size(I));
    goodLbl = [];
    for iLoc = 1:numel(x0)
        lines = round(y0(iLoc))+[-1:1];
        cols = round(x0(iLoc))+[-1:1];
        lines(lines<=0 | lines > size(sproutLocs,1)) = [];
        cols(cols<=0 | cols > size(sproutLocs,2)) = [];
        thisGoodLabel = unique(sproutsLbltmp(lines,cols));
        thisGoodLabel(thisGoodLabel==0) = [];
        goodLbl = unique([goodLbl;thisGoodLabel]);
        sproutLocs(round(y0(iLoc)),round(x0(iLoc))) = 1;
    end
    sproutLocs = logical(imdilate(sproutLocs,strel('disk',5,8))-imdilate(sproutLocs,strel('disk',3,8)));
    
    % remove bad labels
    for iLbl = goodLbl'
        sprouts = sprouts | sproutsLbltmp==iLbl;
    end
    
    sproutsLbl = bwlabel(sprouts);
    sproutsProps = regionprops(sproutsLbl,'pixelList');
    
    for iProps = 1:numel(sproutsProps)
        thesePoints = sproutsProps(iProps).PixelList;
        distMtrx = sqrt((thesePoints(:,1)-thesePoints(:,1)').^2 + (thesePoints(:,2)-thesePoints(:,2)').^2);
        nearestDist = distMtrx(distMtrx <= sqrt(2) & distMtrx~=0);
        sproutsProps(iProps).Length = sum(nearestDist)/2;
    end
    
    % add the image name
    for iProps = 1:numel(sproutsProps)
        sproutsProps(iProps).Image = t;
    end
    
    
    resultTableOneIm = struct2table(sproutsProps);
    % convert from pixels to microns
    resultTableOneIm.Length = resultTableOneIm.Length*imRes;
    resultTableOneIm = removevars(resultTableOneIm,{'PixelList','Image'});
    resultTable = [resultTable; resultTableOneIm];
    
    resultTableSummary(idx).Image = t;
    resultTableSummary(idx).Mean_Length = mean(resultTableOneIm.Length,'omitnan');
    resultTableSummary(idx).Std_Length = std(resultTableOneIm.Length,'omitnan');
    resultTableSummary(idx).n_sprouts = length(resultTableOneIm.Length);
   
    writetable(resultTableOneIm,fullfile(resultFolder,'results.xlsx'),'Sheet',t)
    
    resultImage = uint8(cat(3,Ihisteq+contourDist1+roughContour,Ihisteq+sproutLocs,Ihisteq+imdilate(sprouts,strel('disk',1)))*255);
    sproutImage = uint8(cat(3,Ihisteq+contourDist1+roughContour,Ihisteq+skeleton,Ihisteq+imdilate(sprouts,strel('disk',1)))*255);
    imwrite(resultImage,fullfile(resultFolder,strrep(t,'.tif','-result.tif')))
    imwrite(sproutImage,fullfile(resultFolder,strrep(t,'.tif','-sprouts.tif')))
end


% save data
save(fullfile(resultFolder,'results.mat'),'resultTable')
writetable(struct2table(resultTableSummary),fullfile(resultFolder,'results.xlsx'),'Sheet','Summary')


