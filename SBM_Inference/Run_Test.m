rng(2018)

%% Real Run
path = '~/Documents/hcp_data/HCP_components.mat';
nperms = 10000;
topComps = 809;
dropNets = [-1];
sanity=0;

[sbm.pval,sbm.log_lik_true,sbm.log_lik_null_median] = FundamentalUnits_SBM_Test(path,nperms,topComps,dropNets,sanity);

%% Pre-shuffle the labels and run again as a sanity check
path = '~/Documents/hcp_data/HCP_components.mat';
nperms = 10000;
topComps = 809;
dropNets = [-1];
sanity=1;

[fake.pval,fake.log_lik_true,fake.log_lik_null_median] = FundamentalUnits_SBM_Test(path,nperms,topComps,dropNets,sanity);

save('SBM_Inference.mat','sbm');
