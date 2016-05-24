###################################################
# Section 10.4 Estimating a Table of Means
###################################################

library(LearnBayes)

 data(iowagpa)
 rlabels = c("91-99", "81-90", "71-80", "61-70", "51-60", "41-50", 
     "31-40", "21-30")
 clabels = c("16-18", "19-21", "22-24", "25-27", "28-30")
 gpa = matrix(iowagpa[, 1], nrow = 8, ncol = 5, byrow = T)
 dimnames(gpa) = list(HSR = rlabels, ACTC = clabels)
 gpa

 samplesizes = matrix(iowagpa[, 2], nrow = 8, ncol = 5, byrow = T)
 dimnames(samplesizes) = list(HSR = rlabels, ACTC = clabels)
 samplesizes

 act = seq(17, 29, by = 3)
 matplot(act, t(gpa), type = "l", lwd = 3, 
  xlim = c(17, 34), col=1:8, lty=1:8)
 legend(30, 3, lty = 1:8, lwd = 3, legend = c("HSR=9", "HSR=8", 
     "HSR=7", "HSR=6", "HSR=5", "HSR=4", "HSR=3", "HSR=2"), col=1:8)

S=readline(prompt="Type  <Return>   to continue : ")

 MU = ordergibbs(iowagpa, 5000)

 postmeans = apply(MU, 2, mean)
 postmeans = matrix(postmeans, nrow = 8, ncol = 5)
 postmeans=postmeans[seq(8,1,-1),]
 dimnames(postmeans)=list(HSR=rlabels,ACTC=clabels)
 round(postmeans,2)

windows()
matplot(act, t(postmeans), type = "l", lty=1:8, lwd = 3, col = 1, xlim = c(17, 34))
 legend(30, 3, lty = 1:8, lwd = 2, legend = c("HSR=9", "HSR=8", 
     "HSR=7", "HSR=6", "HSR=5", "HSR=4", "HSR=3", "HSR=2"))

 postsds = apply(MU, 2, sd)
 postsds = matrix(postsds, nrow = 8, ncol = 5)
 postsds=postsds[seq(8,1,-1),]
 dimnames(postsds)=list(HSR=rlabels,ACTC=clabels)
 round(postsds,3)

 s=.65
 se=s/sqrt(samplesizes)
 round(postsds/se,2)

S=readline(prompt="Type  <Return>   to continue : ")

 FIT=hiergibbs(iowagpa,5000)

windows()
 par(mfrow=c(2,1))
 plot(density(FIT$beta[,2]),xlab=expression(beta[2]),
  main="HIGH SCHOOL RANK")
 plot(density(FIT$beta[,3]),xlab=expression(beta[3]),
  main="ACT SCORE")
 quantile(FIT$beta[,2],c(.025,.25,.5,.75,.975))

 quantile(FIT$beta[,3],c(.025,.25,.5,.75,.975))

 quantile(FIT$var,c(.025,.25,.5,.75,.975))

 posterior.means = apply(FIT$mu, 2, mean)
 posterior.means = matrix(posterior.means, nrow = 8, ncol = 5, 
  byrow = T)

S=readline(prompt="Type  <Return>   to continue : ")

windows()
 par(mfrow=c(1,1))
 matplot(act, t(posterior.means), type = "l", lwd = 3, lty=1:8, col=1,
   xlim = c(17, 34))
 legend(30, 3, lty = 1:8, lwd = 2, legend = c("HSR=9", "HSR=8", 
     "HSR=7", "HSR=6", "HSR=5", "HSR=4", "HSR=3", "HSR=2"))

 p=1-pnorm((2.5-FIT$mu)/.65)
 prob.success=apply(p,2,mean)

 prob.success=matrix(prob.success,nrow=8,ncol=5,byrow=T)
 dimnames(prob.success)=list(HSR=rlabels,ACTC=clabels)
 round(prob.success,3)
