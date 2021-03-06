function [options] = defaultOptions_mcbDNMF()
    %[options] = defaultOptions_mcbDNMF()
%     options.thr = 2;
%     options.patchSize = [64 64];
%     options.stride = 56;
%     options.overlapThr = 0.25;
%     options.temporalCorrThr = 0.8;
%     options.minSkew = 1;
%     options.sizeRange = [50 1000];

    options.thr = 0.6;
    options.patchSize = [64 64];
    options.stride = 56;
    options.overlapThr = 0.5;
    options.temporalCorrThr = 0.8;
    options.minSkew = 0;
    options.sizeRange = [30 2000];
            
end