
NumCompCraddock = 100;

fout_mu = zeros(9,1);
fout_std = zeros(9,1);
fout_all_mu = zeros(9,1);
fout_all_std = zeros(9,1);

Nums = [100:100:900];

results = zeros(size(pheno,2),9);

for Num = 1:9
    clear featuremat_craddock
    phenotype_test_predict_craddock = zeros(n,size(pheno,2));
    NUM = num2str(Nums(Num));
    CraddockCorrTemplate = [DataDir 'Connectomes/[Subject]/craddock_' NUM '/corr.mat'];

    featuremat_craddock = load_connectomes(CraddockCorrTemplate,subs);

    if (nFold == 1)
        folds = zeros(n,1);
        folds(Test==1) = 1;
        folds(Train==1) = 2;
        featuremat = featuremat_all;
        pheno = pheno_all;
    end

    train_idx = folds==2;
    test_idx = folds==1;
    
    dims = mledim(featuremat_craddock(train_idx,:)',10,20);
    dims_all(Num) = dims;
    
    [coeff, score, latent, ~, exp] = pca(featuremat_craddock(train_idx,:));

    mu = mean(featuremat_craddock(train_idx,:));
    x = bsxfun(@minus,featuremat_craddock(train_idx,:),mu);
    xt = bsxfun(@minus,featuremat_craddock(test_idx,:),mu);

    A = (pinv(coeff)*x')';
    A_test = (pinv(coeff)*xt')';

    n1t = A_test(:,1:NumCompCraddock) * coeff(:,1:NumCompCraddock)';
    r = fast_corr(xt',n1t');
    fout_mu(Num) = mc_FisherZ(mean(mc_FisherZ(r)),1);
    fout_std(Num) = mc_FisherZ(std(mc_FisherZ(r)),1);
    
    AllComp = 809;
    
    n1t = A_test(:,1:AllComp) * coeff(:,1:AllComp)';
    r = fast_corr(xt',n1t');
    fout_all_mu(Num) = mc_FisherZ(mean(mc_FisherZ(r)),1);
    fout_all_std(Num) = mc_FisherZ(std(mc_FisherZ(r)),1);

    for iPheno = 1:size(pheno,2)
        X = [ones(n_train,1) A(:,1:NumCompCraddock)];
        b = pinv(X'*X)*X'*pheno(train_idx,iPheno);

        phenotype_test_predict_craddock(test_idx,iPheno) = [ones(n_test,1) A_test(:,1:NumCompCraddock)]*b;
    end
    
    temp = [];
    temp = diag(corr(pheno(Test==1,:),phenotype_test_predict_craddock(Test==1,:)));
    results(:,Num) = temp;
end

crad_ROIs = Nums;
crad_dim = dims_all;
crad_pred = results;
crad_fout_mu = fout_mu;
crad_fout_std = fout_std;
crad_fout_all_mu = fout_all_mu;
crad_fout_all_std = fout_all_std;
names = {'genexec','procspeed','pmat','ext','int','attn','O','C','E','A','N'};

save(fullfile(OutputPath,'craddock.mat'),'crad_ROIs','crad_dim','crad_pred','crad_fout_mu','crad_fout_std','crad_fout_all_mu','crad_fout_all_std','names');
