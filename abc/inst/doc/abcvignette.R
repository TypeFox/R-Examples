### R code from vignette source 'abcvignette.Rnw'

###################################################
### code chunk number 1: abcvignette.Rnw:56-57
###################################################
rm(list=ls())


###################################################
### code chunk number 2: abcvignette.Rnw:97-99
###################################################
library(abc)
require(abc.data)


###################################################
### code chunk number 3: abcvignette.Rnw:104-106 (eval = FALSE)
###################################################
## help(package="abc")
## help(abc)


###################################################
### code chunk number 4: abcvignette.Rnw:114-115 (eval = FALSE)
###################################################
## example(abc)


###################################################
### code chunk number 5: abcvignette.Rnw:361-364
###################################################
library(abc)
require(abc.data)
data(musigma2)


###################################################
### code chunk number 6: abcvignette.Rnw:398-401
###################################################
require(abc.data)
data(human)
stat.voight


###################################################
### code chunk number 7: ssplot
###################################################
par(mfcol = c(1,3), mar=c(5,3,4,.5))
boxplot(stat.3pops.sim[,"pi"]~models, main="Mean nucleotide diversity")
boxplot(stat.3pops.sim[,"TajD.m"]~models, main="Mean Tajima's D")
boxplot(stat.3pops.sim[,"TajD.v"]~models, main="Var in Tajima's D")


###################################################
### code chunk number 8: ss
###################################################
par(mfcol = c(1,3), mar=c(5,3,4,.5))
boxplot(stat.3pops.sim[,"pi"]~models, main="Mean nucleotide diversity")
boxplot(stat.3pops.sim[,"TajD.m"]~models, main="Mean Tajima's D")
boxplot(stat.3pops.sim[,"TajD.v"]~models, main="Var in Tajima's D")


###################################################
### code chunk number 9: modsel
###################################################
cv.modsel <- cv4postpr(models, stat.3pops.sim, nval=5, tol=.01, method="mnlogistic")
s <- summary(cv.modsel)


###################################################
### code chunk number 10: cv4postprplot
###################################################
plot(cv.modsel, names.arg=c("Bottleneck", "Constant", "Exponential"))


###################################################
### code chunk number 11: cv4postpr
###################################################
plot(cv.modsel, names.arg=c("Bottleneck", "Constant", "Exponential"))


###################################################
### code chunk number 12: abcvignette.Rnw:501-504
###################################################
mytotal <- length(cv.modsel$cvsamples)/length(unique(models))
myexp <- s$conf.matrix$tol0.01[3,3]
misclasstot <- 1-(sum(s$conf.matrix$tol0.01[1,1],s$conf.matrix$tol0.01[2,2],s$conf.matrix$tol0.01[3,3])/sum(s$conf.matrix$tol0.01))


###################################################
### code chunk number 13: abcvignette.Rnw:526-532
###################################################
modsel.ha <- postpr(stat.voight["hausa",], models, stat.3pops.sim, tol=.05, method="mnlogistic")
modsel.it <- postpr(stat.voight["italian",], models, stat.3pops.sim, tol=.05, method="mnlogistic")
modsel.ch <- postpr(stat.voight["chinese",], models, stat.3pops.sim, tol=.05, method="mnlogistic")
summary(modsel.ha)
summary(modsel.it)
summary(modsel.ch)


###################################################
### code chunk number 14: distanceplot
###################################################
res.gfit.bott=gfit(target=stat.voight["italian",], sumstat=stat.3pops.sim[models=="bott",],
statistic=mean, nb.replicate=100)
plot(res.gfit.bott, main="Histogram under H0")


###################################################
### code chunk number 15: distance
###################################################
res.gfit.bott=gfit(target=stat.voight["italian",], sumstat=stat.3pops.sim[models=="bott",],
statistic=mean, nb.replicate=100)
plot(res.gfit.bott, main="Histogram under H0")


###################################################
### code chunk number 16: abcvignette.Rnw:573-580
###################################################
res.gfit.exp=gfit(target=stat.voight["italian",], sumstat=stat.3pops.sim[models=="exp",],
statistic=mean, nb.replicate=100)
res.gfit.const=gfit(target=stat.voight["italian",], sumstat=stat.3pops.sim[models=="const",],
statistic=mean, nb.replicate=100)
summary(res.gfit.bott)
summary(res.gfit.exp)
summary(res.gfit.const)


###################################################
### code chunk number 17: pcaplot
###################################################
gfitpca(target=stat.voight["italian",], sumstat=stat.3pops.sim, index=models, cprob=.1)


###################################################
### code chunk number 18: pca
###################################################
gfitpca(target=stat.voight["italian",], sumstat=stat.3pops.sim, index=models, cprob=.1)


###################################################
### code chunk number 19: ppcplot
###################################################
require(abc.data)
data(ppc)
mylabels <- c("Mean nucleotide diversity","Mean Tajima's D", "Var Tajima's D")
par(mfrow = c(1,3), mar=c(5,2,4,0))
for (i in c(1:3)){
    hist(post.bott[,i],breaks=40, xlab=mylabels[i], main="")
    abline(v = stat.voight["italian", i], col = 2)
}


###################################################
### code chunk number 20: ppc
###################################################
require(abc.data)
data(ppc)
mylabels <- c("Mean nucleotide diversity","Mean Tajima's D", "Var Tajima's D")
par(mfrow = c(1,3), mar=c(5,2,4,0))
for (i in c(1:3)){
    hist(post.bott[,i],breaks=40, xlab=mylabels[i], main="")
    abline(v = stat.voight["italian", i], col = 2)
}


###################################################
### code chunk number 21: abcvignette.Rnw:674-676
###################################################
stat.italy.sim <- subset(stat.3pops.sim, subset=models=="bott")
head(stat.italy.sim)


###################################################
### code chunk number 22: abcvignette.Rnw:686-687
###################################################
head(par.italy.sim)


###################################################
### code chunk number 23: abcvignette.Rnw:702-708
###################################################
cv.res.rej <- cv4abc(data.frame(Na=par.italy.sim[,"Ne"]), stat.italy.sim, nval=10,
tols=c(.005,.01, 0.05), method="rejection")
cv.res.reg <- cv4abc(data.frame(Na=par.italy.sim[,"Ne"]), stat.italy.sim, nval=10,
tols=c(.005,.01, 0.05), method="loclinear")
summary(cv.res.rej)
summary(cv.res.reg)


###################################################
### code chunk number 24: cv4abcplot
###################################################
par(mfrow=c(1,2), mar=c(5,3,4,.5), cex=.8)
plot(cv.res.rej, caption="Rejection")
plot(cv.res.reg, caption="Local linear regression")


###################################################
### code chunk number 25: cv4abc
###################################################
par(mfrow=c(1,2), mar=c(5,3,4,.5), cex=.8)
plot(cv.res.rej, caption="Rejection")
plot(cv.res.reg, caption="Local linear regression")


###################################################
### code chunk number 26: abcvignette.Rnw:755-758
###################################################
res <- abc(target=stat.voight["italian",], param=data.frame(Na=par.italy.sim[, "Ne"]),
sumstat=stat.italy.sim, tol=0.05, transf=c("log"), method="neuralnet")
res


###################################################
### code chunk number 27: abcvignette.Rnw:769-770
###################################################
summary(res)


###################################################
### code chunk number 28: abchistplot
###################################################
hist(res)


###################################################
### code chunk number 29: abchist
###################################################
hist(res)


###################################################
### code chunk number 30: abcplot
###################################################
par(cex=.8)
plot(res, param=par.italy.sim[, "Ne"])


###################################################
### code chunk number 31: abc
###################################################
par(cex=.8)
plot(res, param=par.italy.sim[, "Ne"])


