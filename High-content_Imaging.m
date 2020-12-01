%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%High-content imaging analysis
%Chunyang Dong
%Tian Lab, UC Davis
%11/22/2019
%
%Summary: This script uses the Canny method for edge detection. Put all 
%your images in one folder. Images must be in .tif file format. Please 
%name images by their well number. 
%
%   Inputs
%       -selected individual images
%       
%
%   Outputs
%       -excel containing mean_intensity, median_intensity, 
%           total_intensity, ROIarea
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
%Opens up a dialogue box to choose files
[file, path] = uigetfile('*.tif', 'Select All Image Files', 'MultiSelect','on');
%Adjust fudge factor so threshold for edge detection is more or less strict
fudgeFactor = 1.5;
conn = 8;
%Index each file and obtain mask of image
for i = 1:numel(file)
    img(:,:,:,i) = imread(fullfile(path,file{i})); 
    [~, threshold] = edge(img(:,:,i), 'Canny');
    BW(:,:,i) = edge(img(:,:,i), 'Canny',threshold*fudgeFactor);
    %Dilate image
    se90 = strel('line',5, 90); 
    se0 = strel('line', 5, 0); 
    BWsdil(:,:,i) = imdilate(BW(:,:,i), [se90 se0]);
    %Remove small objects
    BWC(:,:,i) = bwareaopen(BWsdil(:,:,i), 300);
   %To view the binary mask and corresponding original image use
   %imshowpair(BW(:,:n),img(:,:,n), 'montage') where n is the index of the
   %image you want to view. To view final mask replace BW with BWC. 
   
   %Quantify the intensity values in the masked area
   test = BWC(:,:,i);
   y = test; 
   zed = img(:,:,i);
   gamma = img(:,:,i);
   z = zed; 
   z(y == 0) = 0; 
   % set temp to all values in ROI that aren't zero (not part of the mask)
   temp = z(z~=0);
   %Calculate intensity of regions in the mask: https://www.mathworks.com/help/images/measuring-regions-in-grayscale-images.html
   mean_intensity(:,i) = mean(temp);
   median_intensity(:,i) = median(temp);
   total_intensity(:,i) = sum(temp);
   %Calculate the ROI area
   ROIarea(:,i) = sum(sum(gamma.*z));
    
end

%Export data to excel or txt file 
T = table(file', mean_intensity', median_intensity', total_intensity', ROIarea');
s= path;
%Remove special characters or file will not save properly 
s(regexp(s, '[.,/ ]')) = []
%Saves table to excel file where s is the file name
writetable(T, [s '.xls']); 
