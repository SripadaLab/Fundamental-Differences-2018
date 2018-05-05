function tot_log_lik = SBM_NormLogLik(A,C)
% Given a symmetric adjacency matrix A (with uninteresting diagonal
% values), whose entries are presumed to come from a univariate normal
% distribution, and node membership labels C \in 1, ... K, find the MLE estimators
% of \mu and \sigma^2 for each of the
% K on-diagonal blocks + (K Choose 2) off-diagonal blocks, and then
% return the profile log likelihood of the MLE given the observed A,
% conditional on C
%
% Inputs
%
%       A       n*n adjacency matrix. Entries should be in R and normal-ish
%       C       n-length vector of integers specifying community membership
    tot_log_lik = 0;

    for i=1:max(C)
        for j=i:max(C)
            tot_log_lik = tot_log_lik + calc_log_lik_block(i,j);
        end
    end

    function log_lik = calc_log_lik_block(ci,cj)
    % calculate the log likelihood for block consisting of edges
    % originating with community `ci' and terminating in community
    % `cj' In typical usage, this will be the contribution of the
    % present block to the log likelihood for the entire adjacency
    % matric
        vec = extract_block(ci,cj); % get the edges
        mu = mean(vec);
        s2 = var(vec,1); % MLE of variance

        n = numel(vec);
        log_lik = -n/2 * log(2 * pi) - n/2 *  log(s2) - sum((vec - mu).^2) / 2 * s2;
    end

    function vec = extract_block(ci,cj)
    % extract edges corresponding to block consisting of edges linking
    % community i (ci) to community j (cj)
    %
    % Note: edges will be returned as a vector

    % handle on-diagonal and off-diagonal blocks separately
        if (ci==cj)
            vec = extract_block_on_diag(ci,cj);
        elseif (ci~=cj)
            vec = extract_block_off_diag(ci,cj);
        end
    end

    function vec = extract_block_on_diag(ci,cj)
    % Extract as vector only elements that are above the diagonal.
    % This only makes sense for blocks that are along the diagonal,
    % since they contain redundant entries along with the
    % uninteresting diagonal
        block = A(C==ci,C==cj);
        mask = logical(triu(ones(size(block)),1));
        vec = block(mask);
    end

    function vec = extract_block_off_diag(ci,cj)
    % Extract as vector elements for a given block. This will extract
    % edges originating in community ci and terminating in community
    % cj. If A is symmetric, this will return the same elements
    % (though not in the same order) as calling with cj, ci.
        block = A(C==ci,C==cj);
        vec = block(:);
    end

end
