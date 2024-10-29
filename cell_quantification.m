clc;close all; clear all;warning off;
%% Load and convert
tic;
location = '~/Projects/TNBC/data/*'; %  folder in which your images are
ds = imageDatastore(location);         %  Creates a datastore for all images in your folder
i = 0;
while hasdata(ds) 
    i=i+1;                              % cycle counter
    [rgb,info] = read(ds) ;             % read image from datastore and info
    I = rgb2gray(rgb);                  % convert to gray scale
%% Scale bar measurings
    if i == 1
        %imgTresh = I < 60;
        imgTresh = rgb(:, :, 3) < 90; % filtering only the black of the image
        BW_out = bwpropfilt(imgTresh, 'Area', [600, 1500]); % Taking only the scale bar
        imshow(BW_out)
        stats1 = regionprops(BW_out,'BoundingBox'); % Measuring the scale bar
        pix = stats1.BoundingBox(1,3);  % Take the length
        scale = (200)/pix; % Convert pixels to micormeters
    end
%% Find edges
    [BW,threshold] = edge(I,'sobel');
    fudgeFactor = 0.75;
    BWs = edge(I,'sobel',threshold*fudgeFactor);
%% Dilate
    se90 = strel('line',3,90);
    se0 = strel('line',3,0);
    BWsdil = imdilate(BWs,[se90 se0]);
%% Fill
    BWdfill = imfill(BWsdil,'holes');
    %BWdfill = imfill(BWs,'holes');
%% Remove extra things in borders
    BWnobord = imclearborder(BWdfill,4);
%% Smoothen object
    seD = strel('diamond',1);
    BWfinal = imerode(BWnobord,seD);
    BWfinal = imerode(BWfinal,seD);
    BWfinal = imerode(BWfinal,seD);
    BWfinal = imerode(BWfinal, seD);
    BWfinal = imerode(BWfinal, seD);
    BWfinal = imerode(BWfinal, seD);
%% Extracting name of figure, measuring area and creating circle
    name(i) = convertCharsToStrings(info.Filename(end-11:end-4));
    stats = regionprops('table',BWfinal,'Centroid','Area');
    index = find(stats.Area >= 0.2*max(stats.Area));
    center_pix = stats.Centroid(index,:);
    area_pix = stats.Area(index,:);
    r_pix = (area_pix/pi).^(1/2);
%     figure(2)
    figure(i);
    subplot(1,2,1),imshow(I)
    title(name(i))
    subplot(1,2,2),imshow(imoverlay(I,BWfinal));hold on;
    title(name(i))
    viscircles(center_pix,r_pix);
    number = length(index);
%     r = r_pix*scale;
    area = area_pix*scale^2;
    area_mean = mean(area);
    area_std = std(area);
%     area_median = median(area);
    result(i,:) = [area_mean area_std number];    
end
%% Saving a table with the results in Excel
T = array2table(result,'VariableNames',{'Mean (microm sq)','Standard dev (microm sq)','Number of Cells'},...
    'RowNames',name);
filename = 'results.xlsx';
writetable(T,filename,'WriteRowNames',true,'UseExcel',false)
Time = toc;
% sgtitle(['Total Processing Time was: ',num2str(Time),' seconds'])