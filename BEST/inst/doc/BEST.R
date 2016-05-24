### R code from vignette source 'BEST.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(continue="  ")


###################################################
### code chunk number 2: loadBEST
###################################################
library(BEST)


###################################################
### code chunk number 3: data2grps
###################################################
y1 <- c(5.77, 5.33, 4.59, 4.33, 3.66, 4.48)
y2 <- c(3.88, 3.55, 3.29, 2.59, 2.33, 3.59)


###################################################
### code chunk number 4: data2grpsPriors
###################################################
priors <- list(muM = 6, muSD = 2)


###################################################
### code chunk number 5: run2grps
###################################################
BESTout <- BESTmcmc(y1, y2, priors=priors, parallel=FALSE)


###################################################
### code chunk number 6: meanDiff2grps
###################################################
plot(BESTout)


###################################################
### code chunk number 7: meanDiffGTzero
###################################################
meanDiff <- (BESTout$mu1 - BESTout$mu2)
meanDiffGTzero <- mean(meanDiff > 0)


###################################################
### code chunk number 8: ttest2grps
###################################################
t.test(y1, y2)


###################################################
### code chunk number 9: meanDiff2grpsMore
###################################################
plot(BESTout, compVal=1, ROPE=c(-0.1,0.1))


###################################################
### code chunk number 10: sd2grps
###################################################
plot(BESTout, which="sd")


###################################################
### code chunk number 11: summary2g
###################################################
summary(BESTout)


###################################################
### code chunk number 12: summary2gMore
###################################################
summary(BESTout, credMass=0.8, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
          compValeff=1) 


###################################################
### code chunk number 13: class2g
###################################################
class(BESTout)
print(BESTout)


###################################################
### code chunk number 14: ppd2grps
###################################################
plotPostPred(BESTout)


###################################################
### code chunk number 15: plotAll2grps
###################################################
plotAll(BESTout)


###################################################
### code chunk number 16: attach2grps
###################################################
names(BESTout)
meanDiff <- (BESTout$mu1 - BESTout$mu2)
meanDiffGTzero <- mean(meanDiff > 0)
meanDiffGTzero


###################################################
### code chunk number 17: vars2grps
###################################################
varRatio <- BESTout$sigma1^2 / BESTout$sigma2^2
median(varRatio)
hdi(varRatio)
mean(varRatio > 1)
plotPost(varRatio, xlim=c(0, 30))


###################################################
### code chunk number 18: run1grp
###################################################
y0 <- c(1.89, 1.78, 1.30, 1.74, 1.33, 0.89)
BESTout1g <- BESTmcmc(y0, priors=NULL, parallel=FALSE)


###################################################
### code chunk number 19: mean1grp
###################################################
BESTout1g
plot(BESTout1g)


###################################################
### code chunk number 20: plotAll1grp
###################################################
plotAll(BESTout1g)


###################################################
### code chunk number 21: attach1grp
###################################################
names(BESTout1g)
length(BESTout1g$nu)
variance <- BESTout1g$sigma^2
plotPost(variance, xlim=c(0, 3))


