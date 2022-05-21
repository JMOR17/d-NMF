function [ROIs, Cs, coherence, skew, sz, tElapsed] = mcb_DNMF(path_to_video, options)
    % [ROIs, Cs, coherence, skew, sz, tElapsed] = mcb_DNMF(path_to_video, options)
    fileName = path_to_video;
    Y = bigread2(fileName,1);
%     Y = double(Y);
    
    tStart = tic;
    [szA, szB, nFrames] = size(Y);
    
    if(nargin<2)
        options = defaultOptions_mcbDNMF();
    end
    
    thr = options.thr;
    patchSize = options.patchSize;
    stride = options.stride;
    overlapThr = options.overlapThr;
    temporalCorrThr = options.temporalCorrThr;
    minSkew = options.minSkew;
    sizeRange = options.sizeRange;
    
    [cROIs0, Cs0, coherence0, skew0, sz0] = DNMF_General(Y, thr, patchSize, stride, overlapThr);

    
    valid = skew0>minSkew & isbetween(sz0,sizeRange(1),sizeRange(2));
    cROIs1 = cROIs0(:,valid);
    Cs1 = Cs0(valid,:);
    coherence1 = coherence0(valid);
    skew1 = skew0(valid);
    sz1 = sz0(valid);
    
    ROIs1 = reshape(full(cROIs1),[szA, szB, size(cROIs1,2)]);    

    
    %%
    zCs = zscore(Cs1,[],2);
    cROIs_binary = double(cROIs1>0);
    spatialOverlap_0 = cROIs_binary'*cROIs_binary;
    sizes = sum(cROIs_binary);
    [sz1, sz2] = meshgrid(sizes);
    temp = zeros(size(sz1,1),size(sz1,2),2);
    temp(:,:,1) = sz1;
    temp(:,:,2) = sz2;
    minSize = min(temp,[],3);

    spatialOverlap = spatialOverlap_0./minSize;
    temporalCorr = zCs*zCs'/(size(zCs,2)-1);
    
    linked = (spatialOverlap>0) & (temporalCorr>temporalCorrThr);    
    [SSS,CCC] = graphconncomp(sparse(linked));

    
    ROIs = zeros(szA, szB, max(CCC));
    Cs = zeros(max(CCC), size(Cs1,2));
    coherence = NaN(max(CCC),1);
    skew = NaN(max(CCC),1);
    sz = NaN(max(CCC),1);
    for ii=1:max(CCC)
        these = CCC==ii;        
        ROIs(:,:,ii) = max(ROIs1(:,:,these),[],3);
        Cs(ii,:) = nanmean(Cs1(these,:),1);
        [ch, sk, z] = evaluateROIs(reshape(ROIs(:,:,ii),[],1), Cs(ii,:), [szA szB]);
        coherence(ii) = ch;
        skew(ii) = sk;
        sz(ii) = z;
    end
   
    tElapsed = toc(tStart);
end