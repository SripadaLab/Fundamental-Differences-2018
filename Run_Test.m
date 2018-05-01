path = '~/Documents/hcp_data/HCP_components.mat';
nperms = 100;
topComps = 10;
dropNets = [-1];
sanity=0;

ps_real = FundamentalUnits_SBM_Test(path,nperms,topComps,dropNets,sanity);

path = '~/Documents/hcp_data/HCP_components.mat';
nperms = 100;
topComps = 10;
dropNets = [-1];
sanity=1;

ps_fake = FundamentalUnits_SBM_Test(path,nperms,topComps,dropNets,sanity);
