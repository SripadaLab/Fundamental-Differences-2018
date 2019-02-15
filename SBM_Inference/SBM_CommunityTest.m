function [pval, log_lik_true, log_lik_null_median, log_lik_null] = SBM_CommunityTest(A,C,nperms)
% Conduct permutation testing to assess whether the community assignments in C are ``meaningful''. 
%
% Inputs
%       A       n*n symmetric adjacency matrix, with 1's on the diagonal
%       C       n-length vector of integers specifying community membership
%       nperms  number of shufflings of C to try

    if (~issymmetric(A))
        error('A must be symmetric, but this appears to not be the case')
    end

    log_lik_true = SBM_NormLogLik(A,C);
    log_lik_null = zeros(nperms,1);

    for i=1:nperms
        perm_C = C(randperm(numel(C)));
        log_lik_null(i) = SBM_NormLogLik(A,perm_C);
    end

    pval = sum(log_lik_true <= log_lik_null) / nperms;
    pval = max(pval,1/nperms); % avoid p=0

    log_lik_null_median = median(log_lik_null);
end
