function result = getColonyOutline(img, pixelSize, dZ)
%This function creates the colony cross-sectional outline
rescaleFactor = pixelSize/dZ*3;
result = struct.empty;

% calculate gradient for image
imgrad = zeros(size(img,1),size(img,2),size(img,3));
for k = 1:size(img,3)
    imgrad(:,:,k) = imgradient(img(:,:,k));
end

% invert z-direction of image 
img = img(:,:,end:-1:1);
imgrad = imgrad(:,:,end:-1:1);

% determine threshold
[~,maxSlice] = max(squeeze(sum(sum(img,1),2)));
threshold = multithresh(imgrad(:,:,maxSlice));

% determine mask based on threshold
mask = zeros(size(imgrad,1),size(imgrad,2),3*size(imgrad,3));
mask_increase = zeros(3*size(imgrad,3),1);
volume = 0;
for i = 1:size(imgrad,3)
     mask_slice = imgrad(:,:,i)>threshold;
     mask_slice = imclose(mask_slice,strel('disk',3));
     mask_slice_filled = imfill(mask_slice, 'holes');
     mask_slice_filled = imopen(mask_slice_filled,strel('disk',3));

     mask_increase(3*(i-1)+1:3*(i-1)+3) = sum(mask_slice_filled(:))- sum(mask_slice(:));
     mask(:,:,3*(i-1)+1) = mask_slice_filled;
     mask(:,:,3*(i-1)+2) = mask_slice_filled;
     mask(:,:,3*(i-1)+3) = mask_slice_filled;
     volume = volume + sum(mask_slice_filled(:))*pixelSize/(rescaleFactor)*pixelSize/(rescaleFactor)*dZ;
end

% determine first layer in which the entire colony is filled
firstSlice = find(mask_increase >500,1);

% remove all but largest object
obj = bwconncomp(mask);
stats = regionprops(obj, 'Area');
[~,largestObject] = max([stats.Area]);
mask_new = zeros(size(mask));
mask_new(obj.PixelIdxList{largestObject}) = 1;
mask = mask_new;

% get 2D outline for tilt estimate
mask_2D = sum(mask,3)>0;
[y_ids, x_ids] = ind2sub(size(mask_2D),find(mask_2D));
x_c = mean(x_ids);
y_c = mean(y_ids);

[outline_y, outline_x] = ind2sub(size(mask_2D),find(bwdist(~mask_2D)<2.*mask_2D));
distToCenter = hypot(outline_x - x_c, outline_y - y_c);
angle = asind((outline_x - x_c)./distToCenter);
angleChoice = -89:10:90;
coords = zeros(length(angleChoice),3);
for i = 1:length(angleChoice)
    [~,ind] = min(abs(angle - angleChoice(i)));
    s = squeeze(mask(outline_y(ind), outline_x(ind),:));
    plane = find(s>0,1);
    coords(i,:) = [outline_y(ind), outline_x(ind), plane+1];
end

% determine tilt and reference points
model = pcfitplane(pointCloud(coords),50);
stats = regionprops(mask, 'Centroid');
centroid = stats.Centroid;
centroid = projectPoints([centroid(2) centroid(1) centroid(3)], model);


% determine distance to plane (height) and distance to centroid 

mask(:,:,end+1) = zeros(size(mask,1),size(mask,2));
mask_distToSurface = bwdist(~mask);
mask_outline = mask_distToSurface<2 & mask;

% correct outline for points below first filled slice
mask_dist = zeros(size(mask));
for i = 1:firstSlice
    slice = mask(:,:,i);
    if sum(slice(:))>0
        mask_dist(:,:,i) = bwdist(~slice);
        mask_outline(:,:,i) = zeros(size(slice));
        [outline_y, outline_x] = ind2sub(size(mask_outline),find(mask_dist(:,:,i)<2 & slice));
        distToCenter = hypot(outline_x - x_c, outline_y - y_c);
        angle = round(4*(asind((outline_x - x_c)./distToCenter)+90))/4; 
        for v = 0:0.25:360
            temp = distToCenter;
            temp(angle~=v) = 0;
            [~,ind] = max(temp);
            mask_outline(outline_y(ind), outline_x(ind),i) = 1;
        end

    end
end

[x,y,z] = ind2sub(size(mask),find(mask_outline>0));
[points_projected, height] = projectPoints([x y z], model);
if mean(height)<0
    height = -height;
end
height = pixelSize/rescaleFactor*height;
distanceToCenter = pixelSize/rescaleFactor*sqrt((points_projected(:,1)-centroid(1)).^2+(points_projected(:,2)-centroid(2)).^2+(points_projected(:,3)-centroid(3)).^2);

distBins = 0:1:ceil(max(distanceToCenter));
distVsHeight = zeros(length(distBins-1),2);
for i = 1:length(distBins)-1
    distVsHeight(i,1) = 0.5*(distBins(i)+distBins(i+1));
    indices = distanceToCenter>distBins(i) & distanceToCenter<= distBins(i+1);
    distVsHeight(i,2) = nanmean(height(indices));
    distVsHeight(i,3) = nanstd(height(indices));
    distVsHeight(i,4) = sum(indices);
end

result(end).distVsHeight = distVsHeight;
result(end).volume = volume;



    function [resPoints, lambda] = projectPoints(points, model)
        lengthVec = sqrt(model.Normal(1)^2 + model.Normal(2)^2 + model.Normal(3)^2);
        a = model.Parameters(1)/lengthVec;
        b = model.Parameters(2)/lengthVec;
        c = model.Parameters(3)/lengthVec;
        d = model.Parameters(4)/lengthVec;
        tempVal = points(:,1)*a + points(:,2)*b + points(:,3)*c + d;
        lambda = -tempVal./(a^2 + b^2 + c^2);
        resPoints(:,1) = points(:,1)+lambda*a;
        resPoints(:,2) = points(:,2)+lambda*b;
        resPoints(:,3) = points(:,3)+lambda*c;
    end

end