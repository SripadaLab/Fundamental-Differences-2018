library(ggplot2)
library(R.matlab)

results = R.matlab::readMat('SBM_Inference.mat')

sbm = data.frame(results$sbm[,,1])

sbm$CompNum = 1:nrow(sbm)
sbm$LogLikGap = sbm$log.lik.true - sbm$log.lik.null.median


png('SBM_Likelihood.png',units='in',width=5,height=5,res=300)
ggplot(data=sbm) +
    geom_line(aes(x = CompNum, y=log.lik.true,color='blue'),show.legend=FALSE) +
    geom_line(aes(x = CompNum, y=log.lik.null.median,color='red'),show.legend=FALSE) +
    xlab('Component Number') +
    ylab('Profile Log Likelihood')
dev.off()

png('SBM_Likelihood_Diff.png',units='in',width=5,height=5,res=300)
ggplot(data=sbm) +
    geom_line(aes(x = CompNum, y=LogLikGap,color='yellow'),show.legend=FALSE) +
    xlab('Component Number') +
    ylab('Difference in Profile Log Likelihood')
dev.off()
