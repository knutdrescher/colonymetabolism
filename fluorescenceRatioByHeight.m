
function  result = fluorescenceRatioByHeight(img, img_c2, pixelSize, dZ, invert)
%This function calculates the fluorescence and fluorescent ratio in the
%outermost 30 um of the colony as a function of height
rescaleFactor = 0.2;
dist = 30; %um

result = struct.empty;

% invert images if necessary
if invert
    img = img(:,:,end:-1:1);
    img_c2 = img_c2(:,:,end:-1:1);
end

[mask,~] = getColonyMask_v1_cleanup(imresize(img, rescaleFactor), 1);

relevantSlices = squeeze(sum(sum(mask,1),2))>0;
relevantSlices = find(relevantSlices);
img = img(:,:,relevantSlices);
img_c2 = img_c2(:,:,relevantSlices);


result.distToSubstrate = (1:size(img_c2,3))*dZ;
result.mask = mask>0;

mask = imresize(mask, [size(img_c2,1), size(img_c2,2)]);

% determine diameter per slice
diameter = zeros(size(img_c2,3),1);
for k = 1:length(relevantSlices)
    stats = regionprops(mask(:,:,relevantSlices(k)), 'MajorAxisLength');
    diameter(k) = [stats.MajorAxisLength];
end

result(end).diameters = diameter;

mask_dist_final = false(size(img_c2));

% determine intensity in outer area
intensity_ch1 = zeros(length(relevantSlices),1);
intensity_ch2 = zeros(length(relevantSlices),1);
for k = 1:length(relevantSlices)
    mask_dist_slice = bwdist(~mask(:,:,relevantSlices(k))).*mask(:,:,relevantSlices(k));
    mask_dist_slice = mask_dist_slice < (dist/pixelSize) & mask_dist_slice>0;
    
    mask_dist_final(:,:,k) = mask_dist_slice;
    intensity_ch2(k) = sum(img(mask_dist_final(:,:,k)));
    intensity_ch1(k) = sum(img_c2(mask_dist_final(:,:,k)));
end


result.biomassPerSlice = squeeze(sum(sum(mask,1),2))*pixelSize*pixelSize*dZ;
result.biomassPerSlice_distance_30 = squeeze(sum(sum(mask_dist_final,1),2))*pixelSize*pixelSize*dZ;
result.intensity_ch1 = intensity_ch1;
result.intensity_ch2 = intensity_ch2;
result.intensity_ratio = intensity_ch2./intensity_ch1;

