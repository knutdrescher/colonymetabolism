
function  result = deadCellFractionByHeight(img, img_c2, pixelSize, dZ)
%This function calculates ratio of dead cells in the
%outermost 30 um of the colony as a function of height

result = struct.empty;
rescaleFactor = 0.2;
thresh = 2500;
dist = 30; %?m


img = imresize(img, rescaleFactor);
img_c2 = imresize(img_c2, rescaleFactor);

% invert image
img = img(:,:,end:-1:1);
img_c2 = img_c2(:,:,end:-1:1);

[mask,~] = getColonyMask_v1_cleanup(img, 1);

relevantSlices = squeeze(sum(sum(mask,1),2))>0;
relevantSlices = find(relevantSlices);
img_c2 = img_c2(:,:,relevantSlices);

mask = imresize(mask, [size(img_c2,1), size(img_c2,2)]);

deadCells = img_c2>thresh;
fraction_2 = zeros(size(deadCells,3),1);
for k = 1:length(relevantSlices)
    deadCells(:,:,k) = deadCells(:,:,k).*mask(:,:,relevantSlices(k));
    fraction_2(k) = sum(sum(deadCells(:,:,k)))/sum(sum(mask(:,:,relevantSlices(k))));
end

result(1).distToSubstrate = (1:size(deadCells,3))*dZ;
result(1).fraction_DeadCells_wholeArea = fraction_2;


fraction_2 = zeros(size(deadCells,3),1);
mask_dist_final = zeros(size(deadCells));
nSlices = floor(dist/dZ);
for k = 1:length(relevantSlices)
    mask_dist_slice = bwdist(~mask(:,:,relevantSlices(k))).*mask(:,:,relevantSlices(k));
    mask_dist_slice = mask_dist_slice < (dist/pixelSize) & mask_dist_slice>0;
    if relevantSlices(k)<size(mask,3)
        slices_upper = (relevantSlices(k)+1):min(relevantSlices(k)+nSlices, size(mask,3));
        mask_dist_upper = sum(mask(:,:,slices_upper),3);
        mask_dist_upper = mask(:,:,relevantSlices(k)) & (mask_dist_upper < length(slices_upper));
    else
        mask_dist_upper = zeros(size(mask,1),size(mask,2));
    end
    mask_dist_final(:,:,k) = mask_dist_slice | mask_dist_upper;
    deadCells(:,:,k) = deadCells(:,:,k).*mask_dist_final(:,:,k);
    fraction_2(k) = sum(sum(deadCells(:,:,k)))/sum(sum(mask_dist_final(:,:,k)));
end

result(1).fraction_DeadCells_partialArea = fraction_2;
