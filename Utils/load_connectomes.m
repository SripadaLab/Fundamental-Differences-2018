function featuremat = load_connectomes(CorrTemplate, subs)
    n = size(subs,1);

    %load and z-transform connectivity matrices
    for iSubject = 1:n
        fprintf(1,'Loading subject %d of %d\n',iSubject,n);
        corrpath = strrep(CorrTemplate,'[Subject]',subs(iSubject,:));
        r = load(corrpath,'rMatrix');
        r = r.rMatrix;
        z = mc_FisherZ(r);
        z = mc_flatten_upper_triangle(z);
        if (iSubject==1)
            featuremat = zeros(n,numel(z));
        end
        z(isnan(z)) = 0;
        featuremat(iSubject,:) = z;
    end