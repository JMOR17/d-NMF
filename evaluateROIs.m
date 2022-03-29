function [coherence, skew, roiSize] = evaluateROIs(A, C, patchSize)
    % [coherence, skew] = evaluateROIs(A, C, patchSize)
    
    AA = reshape(full(A), patchSize(1), patchSize(2), size(A,2));
    h = [1 1 1; 1 0 1; 1 1 1];
    coherence = zeros(size(C,1),1);
    for i_a = 1:size(AA,3)
        this = squeeze(AA(:,:,i_a));
        c = corrcoef(this, conv2(this,h,'same'));
        coherence(i_a) = c(1,2);
    end
    skew = skewness(full(C),[],2);

    roiSize = full(sum(A>0))';
end