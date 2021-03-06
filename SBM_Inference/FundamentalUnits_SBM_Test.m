function [ps, real_ll, null_ll_median] = FundamentalUnits_SBM_Test(path,nperms,topComps,dropNets,sanity,seed)
% this script implements the SBM Inference Routine as described in the manuscript

    if(~exist('sanity'))
        sanity=0;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Real Work Starts Here %
    %%%%%%%%%%%%%%%%%%%%%%%%%

    comps = load(path);

    C = comps.nets;
    mask = ~ismember(C,dropNets);
    C = C(mask);

    if (sanity==1) % shuffle the C as a sanity check
        C = C(randperm(numel(C)));
    end

    ps = ones(topComps,1);
    real_ll = zeros(topComps,1);
    null_ll_median = zeros(topComps,1);

    parfor i=1:topComps
        tic;
        fprintf(1,'Currently working on Component %d of %d requested\n',i,topComps);
        rng(seed+i);
        A = comps.components(:,:,i);
        if (istriu(A)) % symmetrize if upper triangular
            A = A + A';
        end
        A = A(mask,mask);
        [ps(i), real_ll(i), null_ll_median(i)] = SBM_CommunityTest(A,C,nperms);
        time = toc;
        fprintf(1,'Permutations for Component %d took %f seconds\n\n',i,time);
    end
end
