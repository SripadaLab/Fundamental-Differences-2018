
%to run on Train/Test split
NumComp = [50 100 150];
nFold = 1;

%to run 10-fold CV in training data 
%NumComp = [5:5:500];
%nFold = 10;

NumComps = size(NumComp,2);

n_train = sum(Train);
n_test = sum(Test);

replace = 1;
if (nFold==n)
    replace = 0;
end
rng('default');
rng(12345);
folds = randsample(1:nFold,n_train,replace);

if (nFold == 1)
    folds = zeros(n,1);
    folds(Test==1) = 1;
    folds(Train==1) = 2;
end

phenotype_test_predict = zeros(n,size(pheno,2),NumComps);

for iFold = 1:nFold
    fprintf(1,'.');
    %find train and test data for this fold
    test_idx = folds==iFold;
    train_idx = ~test_idx;
    
    n_trainf = sum(train_idx);
    n_testf = sum(test_idx);

    %reduce the training data
    [coeff, score, latent, ~, exp] = pca(featuremat(train_idx,:));
    Abig_orig = score;
    icasig = coeff';

    %mean center train, and mean center test with train means
    mu = mean(featuremat(train_idx,:));
    x = bsxfun(@minus,featuremat(train_idx,:),mu);
    xt = bsxfun(@minus,featuremat(test_idx,:),mu);

    %calculate expressions for each subject for train and test
    A = (pinv(icasig')*x')';
    A_test = (pinv(icasig')*xt')';

    for iComp = 1:NumComps
        k = NumComp(iComp);
        
        Abig = A(:,1:k);
        Abig_test = A_test(:,1:k);
        %predicting phenotype
        for iPheno = 1:size(pheno,2)
            X = [ones(n_trainf,1) Abig];
            b = pinv(X'*X)*X'*pheno(train_idx,iPheno);
            
            phenotype_test_predict(test_idx,iPheno,iComp) = [ones(n_testf,1) Abig_test]*b;
        end
        
    end
    
end
fprintf(1,'\n');

%check correlation between actual and predicted phenotypes
fold_corr = zeros(size(pheno,2),size(NumComps,2),nFold);
for i = 1:NumComps
    for iFold = 1:nFold
        fold_corr(:,i,iFold) = diag(corr(pheno(folds==iFold,:),phenotype_test_predict(folds==iFold,:,i)));
    end
end

mean_corr = mc_FisherZ(mean(mc_FisherZ(fold_corr),3),1);
