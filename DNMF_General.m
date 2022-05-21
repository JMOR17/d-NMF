function [cROIs, Cs, coherence, skew, sz] = DNMF_General(V, thr, patchSize, stride, overlapThr)

    height = size(V,1);
    width = size(V,2);

    CORE_OR_RANDOM = 1;
    %% Define ROI Pieces & Activity Traces
            
    % Split video into overlapping spatial patches
    indicesA = 1:stride:(height-patchSize(1))+1;
    indicesB = 1:stride:(width-patchSize(2))+1;

    ROIs = cell(1, 1, length(indicesA)*length(indicesB));
    Cs = cell(length(indicesA)*length(indicesB), 1);
    COHERE = cell(length(indicesA)*length(indicesB), 1);
    SKEW = cell(length(indicesA)*length(indicesB), 1);
    SIZE = cell(length(indicesA)*length(indicesB), 1);

    stamp = zeros(height, width);
    count = 1;
    % For each patch 
    for i_A = indicesA
        aa = i_A:(i_A+patchSize(1)-1);
        for i_B = indicesB
            bb = i_B:(i_B+patchSize(2)-1);
 
            thisV = double(V(aa,bb,:));
            
            if(CORE_OR_RANDOM==1)
                % Find cores from original Y 
                A0 = findCores(dv, thr*max(dv,[],3), [1 1 1],15);
                A0 = reshape(A0,patchSize(1)*patchSize(2),[]);
                a0 = sparse(A0);
                roiAND = a0'*a0;
                roiOR = prod(patchSize) - (1-a0)'*(1-a0);
                roiJAC = roiAND./roiOR;
                [~,groups] = graphconncomp(roiJAC>overlapThr);
                temp = zeros(size(A0,1),max(groups));
                for i_group = 1:max(groups)
                    temp(:,i_group) = sum(A0(:,groups==i_group),2)>0;
                end
                A0 = temp;
            else
                % or Initialize Randomly
                A0 = 50;
            end
            
            defoptions = CNMFSetParms;
            defoptions.d1 = size(thisV,1);
            defoptions.d2 = size(thisV,2);
            defoptions.block_size = patchSize;
            if(isempty(A0))
                continue;
            end
            % Pass through SNMF for A, C estimates
            [A,C,B,F] = sparse_NMF_initialization7(thisV,A0,defoptions);
            if(isempty(A))
                continue;
            end
            % Quality check on coherence of footprint and skew of activity
            % (To be used later)
            [coherence, skew, roiSize] = evaluateROIs(A, C, patchSize);
            
            valid = true(size(skew));
            A = A(:,valid);
            C = C(valid,:);
            coherence = coherence(valid);
            skew = skew(valid);
            roiSize = roiSize(valid);
            
            ccc = zscore(A0)'*zscore(A)/(patchSize(1)*patchSize(2)-1);
            [~,i] = max(ccc,[],2);
            indices = unique(i);
            A = A(:,indices);
            C = C(indices,:);
            coherence = coherence(indices);
            skew = skew(indices);
            roiSize = roiSize(indices);                                    
            
            
            % Future: Re-run on residuals until no new segments are found
            
          
            temp = zeros(height, width, size(A,2));
            AA = reshape(full(A), patchSize(1), patchSize(2), size(A,2));
                                   
            h = [1 1 1; 1 0 1; 1 1 1];
            ch = zeros(size(C,1),1);
            for i_a = 1:size(AA,3)
                this = squeeze(AA(:,:,i_a));
                c = corrcoef(this, conv2(this,h,'same'));
                ch(i_a) = c(1,2);
            end
            temp(aa,bb,:) = AA;        
            ROIs{count} = temp;
            Cs{count} = full(C);
            COHERE{count} = coherence(:);
            SKEW{count} = skew(:);
            SIZE{count} = roiSize(:);
            count = count+1;
            if(~isempty(temp))
                stamp = max(stamp, max(temp,[],3));
                clf;
                subplot(1,2,1);
                imagesc(stamp);
                axis square;
                subplot(2,2,2);
                imagesc(max(zscore(thisV,[],3),[],3));
                axis square;
                subplot(2,2,4);
                imagesc(reshape(sum(A,2),patchSize));
                axis square;
                
                drawnow;
            end
        end
    end
    ROIs = cell2mat(ROIs);
    cROIs = sparse(reshape(ROIs,[],size(ROIs,3)));
    Cs = cell2mat(Cs);
    coherence = cell2mat(COHERE);
    skew = cell2mat(SKEW);
    sz = cell2mat(SIZE);

end