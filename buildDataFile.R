rm(list=ls())

#Set paths for HCP data and other derived data
HCPDir = '/net/parasite/HCP/Behavioral/'
DataDir = '/net/parasite/HCP/Scripts/slab/PCA/FundDiffRepo/Data/'

#Load HCP unrestricted, restricted, and retest files
unrestricted = read.csv(paste0(HCPDir,'unrestricted_mangstad_5_19_2017_13_58_48.csv'))
restricted = read.csv(paste0(HCPDir,'RESTRICTED_mangstad_5_19_2017_13_58_52.csv'))
rt = read.csv(paste0(HCPDir,'unrestricted_mangstad_10_30_2017_8_40_29_retest.csv'))

#Select subjects with complete data
included_subs = unrestricted$Subject[unrestricted$X3T_RS.fMRI_PctCompl==100 & unrestricted$T1_Count>0]

#Select variables from unrestricted, restricted, and retest datasets
goodVarsU = c("Subject","Gender","FS_BrainSeg_Vol","NEOFAC_A","NEOFAC_O","NEOFAC_C","NEOFAC_N","NEOFAC_E","fMRI_3T_ReconVrs","PMAT24_A_CR")
goodVarsR = c("Subject","Age_in_Yrs","Family_ID","ASR_Attn_Pct","ASR_Extn_T","ASR_Intn_T")
goodVarsRT = c("Subject","X3T_RS.fMRI_PctCompl","T1_Count")

unrestricted = subset(unrestricted,select=goodVarsU)
restricted = subset(restricted,select=goodVarsR)
rt = subset(rt,select=goodVarsRT)

rt = rt[rt$X3T_RS.fMRI_PctCompl==100 & rt$T1_Count>0,]
rt$Retest = 1
rt = subset(rt,select=c("Subject","Retest"))

data = merge(unrestricted,restricted,by="Subject")
data = data[data$Subject %in% included_subs,]
data = merge(data,rt,by="Subject",all.x=T)
data$Retest[is.na(data$Retest)] = 0

#Now read in our factor analyzed scores and merge
nih = read.csv(paste0(DataDir,'HCP_phenotypic_factors.csv'))
nih$GenExec = as.numeric(as.character(nih$GenExec))
nih$ProcSpeed = as.numeric(as.character(nih$ProcSpeed))
data = merge(data,nih,by="Subject",all.x=T)

#read summary data on motion and merge
motion = read.csv(paste0(DataDir,'Motion.csv'))
motion = motion[,c("Subject","meanFDpre","SupraThresholdFD")]
data = merge(motion,data,by="Subject")

#Now set up some inclusion checks
data$Include.Motion = data$SupraThresholdFD < 480
data$Include.Behav = !is.na(data$PMAT24_A_CR) & !is.na(data$ProcSpeed) & !is.na(data$GenExec) & !is.na(data$NEOFAC_O) & !is.na(data$NEOFAC_C) & !is.na(data$NEOFAC_E) & !is.na(data$NEOFAC_A) & !is.na(data$NEOFAC_N) & !is.na(data$ASR_Attn_Pct) & !is.na(data$ASR_Extn_T) & !is.na(data$ASR_Intn_T) 
data$Include.Retest = as.numeric(data$Retest == 1)
data$Include.NoRetest = data$Include.Behav & !data$Include.Retest
data$IncludeUnrelated = !duplicated(data$Family_ID)
data = data[data$Include.Motion,]

#Some recoding/etc
data$GenderN = as.numeric(data$Gender)
data$TBVc = data$FS_BrainSeg_Vol ^ (1/3)
data$o = 1
data$fMRI_3T_ReconVrs_num = as.numeric(data$fMRI_3T_ReconVrs)
data$Family_ID = as.integer(as.factor(data$Family_ID))


#Now look at families in order to ensure 100 unrelated subjects in test set
length(unique(data$Family_ID))
t = table(data$Family_ID)
d = data.frame(fam=rownames(t),children=c(t))
aggregate(fam~children, d, function(x) length(unique(x)))

temp = data[data$Include.NoRetest,]
length(unique(temp$Family_ID))

t = table(temp$Family_ID)
d = data.frame(fam=rownames(t),children=c(t))
a = aggregate(fam~children,d,function(x) length(unique(x)))

fams = d$fam[d$children==1]
set.seed(12345)
samp = fams[sample(length(fams),100,replace=FALSE)]

data$Include.Test = as.numeric(data$Family_ID %in% samp & data$Include.NoRetest & data$Include.Behav)
data$Include.Train = as.numeric(data$Include.NoRetest & !data$Include.Test)

#check Ns
sum(data$Include.Retest)
sum(data$Include.Test)
sum(data$Include.Train)

#double check to confirm no overlapping subjects
table(data$Include.Retest,data$Include.Test,data$Include.Train)

data1 = subset(data,subset=data$Include.Train | data$Include.Test)

write.csv(data1,paste0(DataDir,'phenotype_910.csv'),row.names=FALSE)

data2 = subset(data,subset=data$Include.Retest==1)
write.csv(data2,paste0(DataDir,'phenotype_rt38.csv'),row.names=FALSE)
