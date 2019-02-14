AnalysisDir = fileparts(mfilename('fullpath'));
%AnalysisDir = pwd
addpath([AnalysisDir '/Utils/']);

DataDir = '/net/parasite/HCP/derivatives/FundDiffRepo/';

PhenotypeFile = [AnalysisDir '/Data/phenotype_910.csv'];

Nodes = 264;
Edges = (Nodes*(Nodes-1))/2;

NumComp = 75;
ControlAll = 0;
DMNTPN = 1;
NumPerms = 10000;

OutputPath = [AnalysisDir '/Results/'];

Model = 's6_Power264_p50f0b_nonaggr_p35mask';
ParamTemplate = [DataDir 'Connectomes/100206/power_264/parameters.mat'];
CorrTemplate = [DataDir 'Connectomes/[Subject]/power_264/corr.mat'];

dat = readtable(PhenotypeFile);
subs = num2str(dat.Subject);

FD = dat.meanFDpre;
gender = dat.GenderN;
TBV = dat.TBVc;
Age = dat.Age_in_Yrs;
Recon = dat.fMRI_3T_ReconVrs_num;

pheno = [dat.GenExec dat.ProcSpeed dat.PMAT24_A_CR dat.ASR_Extn_T ...
         dat.ASR_Intn_T dat.ASR_Attn_Pct dat.NEOFAC_O dat.NEOFAC_C dat.NEOFAC_E ...
         dat.NEOFAC_A dat.NEOFAC_N];

Train = dat.Include_Train;
Test = dat.Include_Test;

n = size(subs,1);
p = (Nodes*(Nodes-1))/2;

featuremat = load_connectomes(CorrTemplate,subs,Nodes);

if (ControlAll)
    nuisance = [Age Age.^2 meanFD meanFD.^2 gender TBV TBV.^2 Recon==2 Recon==3];
    N = size(nuisance,1);
    refnuisance = [mean(Age)*ones(N,1) (mean(Age).^2)*ones(N,1) zeros(N,3) ...
                   mean(TBV)*ones(N,1) (mean(TBV).^2)*ones(N,1) zeros(N,2)];
else
    nuisance = 0*Age;
    refnuisance = nuisance;
end

%need to load retest and _rt data also
RetestPhenotypeFile = [AnalysisDir '/Data/phenotype_rt38.csv'];
RetestCorrTemplate = [DataDir 'Connectomes/[Subject]_rt/power_264/corr.mat'];

dat_rt = readtable(RetestPhenotypeFile);
subs_rt = num2str(dat_rt.Subject);

featuremat_rt1 = load_connectomes(CorrTemplate,subs_rt,Nodes);
featuremat_rt2 = load_connectomes(RetestCorrTemplate,subs_rt,Nodes);
