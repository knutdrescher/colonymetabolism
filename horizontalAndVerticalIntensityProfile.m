function result = horizontalAndVerticalIntensityProfile(img, pixelSize, dZ)
%This function determines the horizontal and vertical fluorescent intensity
%profile as a function of distance to colony surface

rescaleFactor = 0.2;

result = struct.empty;


% resize image for faster processing
img = imresize(img, rescaleFactor);

% calculate slice-wise gradient
imgrad = zeros(size(img));
for i = 1:length(info)
    im_slice = imread(img_path, i);
    im_slice = imresize(im_slice, rescaleFactor);
    img(:,:,i) = im_slice;
    imgrad(:,:,i) = imgradient(im_slice);
end

[~,maxSlice] = max(squeeze(sum(sum(img,1),2)));
threshold = multithresh(imgrad(:,:,maxSlice));

% determine mask based on threshold
mask = zeros(size(imgrad,1),size(imgrad,2),size(imgrad,3));
for i = 1:size(imgrad,3)
    mask_slice = imgrad(:,:,i)>threshold;
    mask_slice_filled = imfill(mask_slice, 'holes');
    if sum(mask_slice_filled(:))> sum(mask_slice(:))+50
         mask_slice_filled = bwareafilt(mask_slice_filled,1);
         mask(:,:,i) = mask_slice_filled;
    end
end

mask_area = squeeze(sum(sum(mask,1),2));
firstSlice = find(mask_area >500,1);
topSlice = find(mask_area >100,1,'last');

[~,largest] = max(mask_area);

% 3D version
mask(:,:,1:firstSlice-1) = 1;

% get central area
mask_2D = sum(mask,3)>0;
[x,y] = ind2sub( size(mask_2D),find(mask_2D));
x_central = mean(x);
y_central = mean(y);
mask_2D = false(size(mask_2D));
mask_2D(round(x_central), round(y_central)) = 1;
mask_2D = imdilate(mask_2D, strel('disk',10));


intensityDistribution = zeros(topSlice-firstSlice+1,4);
for i = firstSlice:topSlice
   intensities = img(:,:,i);
   intensities = intensities(mask_2D);
   intensityDistribution(i-firstSlice+1,1) = (i-firstSlice)*dZ;
   intensityDistribution(i-firstSlice+1,2) = mean(intensities);
   intensityDistribution(i-firstSlice+1,3) = std(intensities);
   intensityDistribution(i-firstSlice+1,4) = numel(intensities);
end

result.intensityDistribution_vertical = intensityDistribution;

% horizontal values - rings from outside to inside
distance_mat = bwdist(~mask(:,:,largest+2));

step = 3;
dist_vals = 0:step:ceil(max(distance_mat(:)));
intensityDistribution = zeros(topSlice-firstSlice+1,4);
for i = 2:length(dist_vals)
    intensities = img(:,:,largest+2);
    sub = dist_vals(i-1)<=distance_mat & distance_mat<dist_vals(i);
    intensities = intensities(sub);
   intensityDistribution(i,1) = (dist_vals(i-1)+step/2)*pixelSize/rescaleFactor;
   intensityDistribution(i,2) = mean(intensities);
   intensityDistribution(i,3) = std(intensities);
   intensityDistribution(i,4) = numel(intensities);
end

result.intensityDistribution_horizontal = intensityDistribution;
result.mask = mask;
        

