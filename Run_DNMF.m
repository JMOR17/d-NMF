clear;
addpath(genpath('.'));
%% Set file path
% files = {   '.\Samples\Sample1\MotionCorrected.tif';
%             '.\Samples2\Sample2\MotionCorrected.tif';           
%             };

% files = {'F:\DNMF\Data\885\885D5._Tsub_mean.tif'};
files = {'F:\Data\MotionCorrected\885\885D5_MC.tif'};
%% Set options
options.thr = 0.6;                    % Threshold for active pixels
options.patchSize = [64 64];        % Size of image patches
options.stride = 56;                % Offset of image patches
options.overlapThr = 0.5;          % Spatial overlap merge threshold 
options.temporalCorrThr = 0.8;      % Temporal correlation merge threshold 
options.minSkew = 0;                % Minimum skew of temporal trace of valid ROIs
options.sizeRange = [30 2000];      % Allowable size range of valid ROIs
        
for i_file = 1:length(files)
    thisFile = files{i_file};
    
    [ROIs, Cs, coherence, skew, sz, tElapsed] = mcb_DNMF(thisFile, options);
    
    folder = fileparts(thisFile);
    outFolder = folder;
    if(~exist(outFolder,'dir'))
        mkdir(outFolder);
    end
    save(fullfile(outFolder,'DNMF_Out.mat'), 'ROIs', 'Cs', 'coherence', 'skew', 'sz', 'options', 'tElapsed', '-v7.3');
        
    outputFolder = fullfile(outFolder,'ROI_Set\');
    if(~exist(outputFolder,'dir'))
        mkdir(outputFolder);
    end
    [~,max_idx] = max(Cs,[],2);
    for i_roi = 1:size(ROIs,3)        
        [a,b] = find(imdilate(ROIs(:,:,i_roi)>0,ones(3)));
        c = boundary(a,b,0.95);

        writeImageJROI_3([a(c) b(c)], 4, max_idx(i_roi), sprintf('r%04d',i_roi), outputFolder);
    end

    zip(strrep(outputFolder,'ROI_Set\','RoiSet_Auto'),'*',outputFolder);
    rmdir(outputFolder, 's');
end