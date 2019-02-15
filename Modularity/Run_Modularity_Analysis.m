path = '~/Documents/hcp_data/HCP_components.mat';
topComps = 100;
dropNets = [-1];

[Q_Fixed, Q_Optimal, C_Optimal] = FundamentalUnits_Modularity_Analysis(path,topComps,dropNets);
