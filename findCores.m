function [ROI] = findCores(V, thr)

    if(nargin<2)
        thr = 2;
    end

    height = size(V,1);
    width = size(V,2);

    CC = bwconncomp(medfilt3(V>thr));
    sz = cellfun(@(x) length(x), CC.PixelIdxList);
    segs = CC.PixelIdxList(sz>50);

    ROIs = zeros(height,width,length(segs));
    for jj=1:length(segs)
        temp = mod(segs{jj}-1,height*width)+1;
        temp2 = zeros(height,width);
        temp2(temp) = 1;
        ROIs(:,:,jj) = temp2;    
    end
    sz = squeeze(sum(sum(ROIs,1),2));
    valid = sz>25;
    ROI = ROIs(:,:,valid);
    
end