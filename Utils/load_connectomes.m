function featuremat = load_connectomes(CorrTemplate, subs, Nodes)
    n = size(subs,1);
    p = (Nodes*(Nodes-1))/2;
    
    featuremat = zeros(n,p);

    %load and z-transform connectivity matrices
    for iSubject = 1:n
        fprintf(1,'Loading subject %d of %d\n',iSubject,n);
        corrpath = strrep(CorrTemplate,'[Subject]',subs(iSubject,:));
        r = load(corrpath,'rMatrix');
        r = r.rMatrix;
        z = mc_FisherZ(r);
        z = mc_flatten_upper_triangle(z);
        featuremat(iSubject,:) = z;
    end