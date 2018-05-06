function [Q_Fixed, Q_Optimal, C_Optimal] = FundamentalUnits_Modularity_Analysis(path,topComps,dropNets)
% this script loops over components and
%       1 - Calculates modularity for the a priori community assignment scheme
%       2 - Runs modularity_und from BCT to find optimal Q and C
%
% NOTE: This algorithm will remove any nodes corresponding to
% dropNets PRIOR to performing all analyses. As a result, Q_Fixed
% and Q_Optimal values correspond to the resultant subgraph, and
% the dimension of C_Optimal may be smaller than the number of ROIs

    comps = load(path);

    C = comps.nets;
    mask = ~ismember(C,dropNets);
    C = C(mask);
    nROI = numel(C); % count VALID rois

    Q_Fixed = zeros(topComps,1);
    Q_Optimal = zeros(topComps,1);
    C_Optimal = zeros(topComps,nROI);

    for i=1:topComps
        tic;
        fprintf(1,'Currently working on Component %d of %d requested\n',i,topComps);
        A = comps.components(:,:,i);
        if (istriu(A)) % symmetrize if upper triangular
            A = A + A';
        end
        A = A(mask,mask);
        Ap = max(0,A); % positive version of A
        An = max(0,-A); % negative version of A (sign-flipped)
        Qp = calc_modularity_und(Ap,C);
        Qn = calc_modularity_und(An,C);
        Q_Fixed(i) = max(Qp,Qn);
        %        Q_Fixed(i) = calc_modularity_und(A,C);
        [~, Q_Opt_p] = modularity_und(Ap);
        [~, Q_Opt_n] = modularity_und(An);
        Q_Optimal(i) = max(Q_Opt_p,Q_Opt_n);
        time = toc;
        fprintf(1,'Modularity Analysis for Component %d took %f seconds\n\n',i,time);
    end
end
