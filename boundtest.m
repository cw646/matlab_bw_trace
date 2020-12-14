
%% Variables
clear all;
folder = cd;
startImage = 1;
numImages = 31;
pixelLength = 20/372;
fps = 3;
timeIncrement = 1/fps;
time = timeIncrement*startImage;

%% Read Files
for i = startImage:numImages
    pngFilename  = sprintf('Layer %d.png', i);
    
    I    = imread(pngFilename);
    Ismooth = imgaussfilt(I,5); % Gaussain smoothing of image
    BW = im2bw(Ismooth); % convert to binary matrix
    B = bwareaopen(BW, 1); % Reduce noise in image
    
    imshow(B)
    
    %imshow(B)
    [B1] = bwboundaries(B);
    
    for j= 1:length(B1)
        bdary = B1{j};
        numPixels(j) = size(bdary, 1);
    end
    
    maxpixel = max(numPixels);
    
    for k= 1:length(B1)
        if numPixels(k) == maxpixel
            sizebound(i) = maxpixel;
            flowboundary = B1{k};
        end
        
    end
    
    
    outputFilename = sprintf('btraceOUT%d.txt', i);
    
    fileID = fopen(outputFilename,'w');
    fprintf(fileID,'%f,%f\n', [flowboundary(:,1),flowboundary(:,2)]');
    fclose(fileID);
    
end
%% Create Boundaries
for  i = startImage:numImages
    infileName = sprintf('btraceOUT%d.txt', i);
    infile = infileName;
    delimiterIn = ',';
    
    A{i} = importdata(infile,delimiterIn);
    
   % peaks ----------------------------------------------
    [imx,imy] = size(BW); 
    
    % tempA{i} =  imy - A{i}(:,2); %for invert
    tempA{i} =  A{i}(:,2);
    
    [peaks{i},locs{i}] = findpeaks(tempA{i}(:,1),1);
    j = length(peaks{i});
    
    for k = 1:j
        peaks{i}(k,2) = A{i}(locs{i}(k,1),1);
        peaks{i}(k,3) = A{i}(locs{i}(k,1),2);
    end  
    
   %------------------------------------------------------
   
   
    A{i}(:,3) = time;
    
    time = time + timeIncrement; %update time
    
end
%% Output To array
L = A{startImage};
for i = startImage:numImages
   % hold on
    %   plot(A{i}(:,2),A{i}(:,1),A{i}(:,3));
    %   plot(A{i}(:,2),A{i}(:,1),'.')
    L = [L;A{i}];
end

%% Mesh grid function

[xi,yi] = meshgrid(1:imy,1:imx);
zi = griddata(L(:,2),L(:,1),L(:,3),xi,yi);


%Select area to exclude data outside of
xb = L(:,2);
yb = L(:,1);
btemp = boundary(xb,yb,1);

in = inpolygon(xi,yi,xb(btemp),yb(btemp)); % identifies the points within the polygon
out = ~in;
zi(out) = NaN;

[Dx,Dy] = gradient(zi);


%   3D Surface:
%   surf(xi,yi,zi)
%   shading interp
%   colormap('winter')
%   colorbar

%% Figure

%Time scale image
im = imagesc(abs(zi));
set(im,'AlphaData',~isnan(zi))
bar = colorbar;
colormap(flipud(parula))
ylabel(bar, 'Time(Seconds)')

%Contour ------------% SET CONTOURS
 hold all
 [C,o] = contour(zi, [1 2 3], 'linewidth', 0.1, ...
     'linecolor','k','ShowText','on');

plot(xb(btemp),yb(btemp),'k','LineWidth',1);
hold on
fill(A{startImage}(:,2), A{startImage}(:,1),'y')

% for i = startImage:numImages
%     hold on
%     plot(peaks{i}(:,3),peaks{i}(:,2),'r*')
% end

% hold o
%q = quiver(xi,yi,-Dx,-Dy);
%q.AutoScaleFactor = 20;
% for i = startLayer:numImages
%     hold on
%     if i > startLayer
%         b = quiver(dy{i}(:,4), dy{i}(:,3),dy{i}(:,1),dy{i}(:,2));
%
%     end
% end