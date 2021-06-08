function [mask,img] = getColonyMask(img, filled)
%This function determines the colony masks based on fluorescence
%Variable filled determines if the output result contains only the
%fluorescent volume or the whole colony

        imgrad = img;

        for i = 1:size(img,3)
            imgrad(:,:,i) = imgradient(img(:,:,i));
        end

        [~,maxSlice] = max(squeeze(sum(sum(img,1),2)));
        % In some cases, this threshold needs to be adapted depending on
        % signal to noise ratio and other imaging conditions
        threshold = 1.5*multithresh(imgrad(:,:,maxSlice));
   
        mask = zeros(size(imgrad,1),size(imgrad,2),size(imgrad,3));
        mask_increase = zeros(size(imgrad,3),1);
        for i = 1:size(imgrad,3)
            mask_slice = imgrad(:,:,i)>threshold;

            mask_slice_open = imopen(mask_slice,strel('disk',2));
            if sum(mask_slice_open(:))==0
                mask_slice = zeros(size(mask_slice));
            end
            mask_slice = imclose(mask_slice,strel('disk',5));
            mask_slice_filled = mask_slice;

            mask_slice_filled = bwareaopen(mask_slice_filled,1000);
            mask(:,:,i) = mask_slice_filled;
            mask_increase(i) = sum(mask_slice_filled(:))- sum(mask_slice(:));

            if filled
                mask(:,:,i) = bwconvhull(mask(:,:,i));
            end
        end
        
        firstSlice = find(mask_increase >500,1);

        mask(:,:,1:firstSlice-1) = 0;
        
        mask = bwareaopen(mask, 50000);


end

