### R code from vignette source 'eRm.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: eRm.Rnw:634-637
###################################################
library("eRm")
res.rasch <- RM(raschdat1)
pres.rasch <- person.parameter(res.rasch)


###################################################
### code chunk number 2: eRm.Rnw:640-642
###################################################
lrres.rasch <- LRtest(res.rasch, splitcr = "mean")
lrres.rasch


###################################################
### code chunk number 3: plotGOF-lrres-rasch (eval = FALSE)
###################################################
## plotGOF(lrres.rasch, beta.subset=c(14,5,18,7,1), tlab="item", conf=list(ia=FALSE,col="blue",lty="dotted"))


###################################################
### code chunk number 4: plotGOF-lrres-rasch-plot
###################################################
plotGOF(lrres.rasch, beta.subset=c(14,5,18,7,1), tlab="item", conf=list(ia=FALSE,col="blue",lty="dotted"))


###################################################
### code chunk number 5: eRm.Rnw:665-668
###################################################
W <- matrix(c(1,2,1,3,2,2,2,1,1,1),ncol=2)
res.lltm <- LLTM(lltmdat2, W)
summary(res.lltm)


###################################################
### code chunk number 6: eRm.Rnw:681-684
###################################################
data(pcmdat2)
res.rsm <- RSM(pcmdat2)
thresholds(res.rsm)


###################################################
### code chunk number 7: plotICC-res-rsm (eval = FALSE)
###################################################
## plotICC(res.rsm, mplot=TRUE, legpos=FALSE,ask=FALSE)


###################################################
### code chunk number 8: plotICC-res-rsm-plot
###################################################
plotICC(res.rsm, mplot=TRUE, legpos=FALSE,ask=FALSE)


###################################################
### code chunk number 9: plotPImap-res-pcm (eval = FALSE)
###################################################
## res.pcm <- PCM(pcmdat2)
## plotPImap(res.pcm, sorted = TRUE)


###################################################
### code chunk number 10: plotPImap-res-pcm-plot
###################################################
res.pcm <- PCM(pcmdat2)
plotPImap(res.pcm, sorted = TRUE)


###################################################
### code chunk number 11: eRm.Rnw:714-716
###################################################
pres.pcm <- person.parameter(res.pcm)
itemfit(pres.pcm)


###################################################
### code chunk number 12: eRm.Rnw:720-724
###################################################
lr<- 2*(res.pcm$loglik-res.rsm$loglik)
df<- res.pcm$npar-res.rsm$npar
pvalue<-1-pchisq(lr,df)
cat("LR statistic: ", lr, "  df =",df, "  p =",pvalue, "\n")


###################################################
### code chunk number 13: eRm.Rnw:741-742
###################################################
grouplpcm <- rep(1:2, each = 10)


###################################################
### code chunk number 14: eRm.Rnw:746-748
###################################################
reslpcm <- LPCM(lpcmdat, mpoints = 2, groupvec = grouplpcm, sum0 = FALSE)
model.matrix(reslpcm)


###################################################
### code chunk number 15: eRm.Rnw:751-752
###################################################
coef(reslpcm, parm="eta")


