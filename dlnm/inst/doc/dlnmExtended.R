### R code from vignette source 'dlnmExtended.Rnw'

###################################################
### code chunk number 1: dlnmExtended.Rnw:51-52
###################################################
options(continue="  ")


###################################################
### code chunk number 2: datadrug
###################################################
library(dlnm)
head(drug, 3)


###################################################
### code chunk number 3: datanested
###################################################
head(nested, 4)


###################################################
### code chunk number 4: Qdrug
###################################################
Qdrug <- as.matrix(drug[,rep(7:4, each=7)])
colnames(Qdrug) <- paste("lag", 0:27, sep="")
Qdrug[1:3,1:14]


###################################################
### code chunk number 5: Qnest
###################################################
Qnest <- t(apply(nested, 1, function(sub) exphist(rep(c(0,0,0,sub[5:14]), 
  each=5), sub["age"], lag=c(3,40))))
colnames(Qnest) <- paste("lag", 3:40, sep="")
Qnest[1:3,1:11]


###################################################
### code chunk number 6: cbdrug
###################################################
cbdrug <- crossbasis(Qdrug, lag=27, argvar=list("lin"),
  arglag=list(fun="ns",knots=c(9,18)))


###################################################
### code chunk number 7: summarydrug
###################################################
summary(cbdrug)


###################################################
### code chunk number 8: pdrug
###################################################
mdrug <- lm(out~cbdrug+sex, drug)
pdrug <- crosspred(cbdrug, mdrug, at=0:20*5)


###################################################
### code chunk number 9: pdrugest1
###################################################
with(pdrug,cbind(allfit,alllow,allhigh)["50",])


###################################################
### code chunk number 10: pdrugest2
###################################################
pdrug$matfit["20","lag3"]


###################################################
### code chunk number 11: plotdrugnoeval (eval = FALSE)
###################################################
## plot(pdrug, zlab="Effect", xlab="Dose", ylab="Lag (days)")
## plot(pdrug, var=60, ylab="Effect at dose 60", xlab="Lag (days)", ylim=c(-1,5))
## plot(pdrug, lag=10, ylab="Effect at lag 10", xlab="Dose", ylim=c(-1,5))


###################################################
### code chunk number 12: plotdrug3d
###################################################
plot(pdrug, zlab="Effect", xlab="Dose", ylab="Lag (days)")


###################################################
### code chunk number 13: plotdruglag
###################################################
plot(pdrug, var=60, ylab="Effect at dose 60", xlab="Lag (days)", ylim=c(-1,5))


###################################################
### code chunk number 14: plotdrugvar
###################################################
plot(pdrug, lag=10, ylab="Effect at lag 10", xlab="Dose",ylim=c(-1,5))


###################################################
### code chunk number 15: cbnested
###################################################
cbnest <- crossbasis(Qnest, lag=c(3,40), argvar=list("bs",degree=2,df=3),
  arglag=list(fun="ns",knots=c(10,30),intercept=F))


###################################################
### code chunk number 16: pnest
###################################################
library(survival)
mnest <- clogit(case~cbnest+strata(riskset), nested)
pnest <- crosspred(cbnest,mnest, cen=0, at=0:20*5)


###################################################
### code chunk number 17: plotnestnoeval (eval = FALSE)
###################################################
## plot(pnest, zlab="OR", xlab="Exposure", ylab="Lag (years)")
## plot(pnest, var=50, ylab="OR for exposure 50", xlab="Lag (years)", xlim=c(0,40))
## plot(pnest, lag=5, ylab="OR at lag 5", xlab="Exposure", ylim=c(0.95,1.15))


###################################################
### code chunk number 18: plotnest3d
###################################################
plot(pnest, zlab="OR", xlab="Exposure", ylab="Lag (years)")


###################################################
### code chunk number 19: plotnestlag
###################################################
plot(pnest, var=50, ylab="OR for exposure 50", xlab="Lag (years)", xlim=c(0,40))


###################################################
### code chunk number 20: plotnestvar
###################################################
plot(pnest, lag=5, ylab="OR at lag 5", xlab="Exposure", ylim=c(0.95,1.15))


###################################################
### code chunk number 21: mylog
###################################################
mylog <- function(x) log(x+1)


###################################################
### code chunk number 22: cbnest2
###################################################
cbnest2 <- crossbasis(Qnest, lag=c(3,40), argvar=list("mylog"),
  arglag=list(fun="ns",knots=c(10,30),intercept=F))
summary(cbnest2)


###################################################
### code chunk number 23: pnest2
###################################################
mnest2 <- clogit(case~cbnest2+strata(riskset), nested)
pnest2 <- crosspred(cbnest2, mnest2, cen=0, at=0:20*5)
plot(pnest2, zlab="OR", xlab="Exposure", ylab="Lag (years)")
plot(pnest2, var=50, ylab="OR for exposure 50", xlab="Lag (years)", xlim=c(0,40))
lines(pnest, var=50, lty=2)
plot(pnest2, lag=5, ylab="OR at lag 5", xlab="Exposure", ylim=c(0.95,1.15))
lines(pnest, lag=5, lty=2)


###################################################
### code chunk number 24: plotnest3d2
###################################################
plot(pnest2, zlab="OR", xlab="Exposure", ylab="Lag (years)")


###################################################
### code chunk number 25: plotnestlag2
###################################################
plot(pnest2, var=50, ylab="OR for exposure 50", xlab="Lag (years)", xlim=c(0,40))
lines(pnest, var=50, lty=2)


###################################################
### code chunk number 26: plotnestvar2
###################################################
plot(pnest2, lag=5, ylab="OR at lag 5", xlab="Exposure", ylim=c(0.95,1.15))
lines(pnest, lag=5, lty=2)


###################################################
### code chunk number 27: fdecay
###################################################
fdecay <- function(x,scale=5) {
  basis <- exp(-x/scale)
  attributes(basis)$scale <- scale
  return(basis)
}


###################################################
### code chunk number 28: cbdrug2
###################################################
cbdrug2 <- crossbasis(Qdrug, lag=27, argvar=list("lin"),
  arglag=list(fun="fdecay",scale=6))
summary(cbdrug2)


###################################################
### code chunk number 29: pdrug2
###################################################
mdrug2 <- lm(out~cbdrug2+sex, drug)
pdrug2 <- crosspred(cbdrug2, mdrug2, at=0:20*5)
plot(pdrug2, zlab="Effect", xlab="Dose", ylab="Lag (days)")
plot(pdrug2, var=60, ylab="Effect at dose 60", xlab="Lag (days)", ylim=c(-1,5))
lines(pdrug, var=60, lty=2)
plot(pdrug2, lag=10, ylab="Effect at lag 10", xlab="Dose", ylim=c(-1,5))
lines(pdrug, lag=10, lty=2)


###################################################
### code chunk number 30: plotdrug3d2
###################################################
plot(pdrug2, zlab="Effect", xlab="Dose", ylab="Lag (days)")


###################################################
### code chunk number 31: plotdruglag2
###################################################
plot(pdrug2, var=60, ylab="Effect at dose 60", xlab="Lag (days)", ylim=c(-1,5))
lines(pdrug, var=60, lty=2)


###################################################
### code chunk number 32: plotdrugvar2
###################################################
plot(pdrug2, lag=10, ylab="Effect at lag 10", xlab="Dose",ylim=c(-1,5))
lines(pdrug, lag=10, lty=2)


###################################################
### code chunk number 33: hist
###################################################
expnested <- rep(c(10,0,13), c(5,5,10))
hist <- exphist(expnested, time=length(expnested), lag=c(3,40))
hist


###################################################
### code chunk number 34: predhist
###################################################
pnesthist <- crosspred(cbnest2, mnest2, at=hist)
with(pnesthist, c(allRRfit,allRRlow,allRRhigh))


###################################################
### code chunk number 35: exp
###################################################
expdrug <- rep(c(10,50,0,20),c(2,1,1,2)*7)


###################################################
### code chunk number 36: fexphist2
###################################################
dynhist <- exphist(expdrug, lag=27)


###################################################
### code chunk number 37: pdynnest
###################################################
pdyndrug <- crosspred(cbdrug2, mdrug2, at=dynhist)


###################################################
### code chunk number 38: plotdynnest
###################################################
plot(pdyndrug,"overall", ylab="Effect", xlab="Time (days)", ylim=c(-10,27), 
  xlim=c(1,50), yaxt="n")
axis(2, at=-1:5*5)
par(new=TRUE)
plot(expdrug, type="h", xlim=c(1,50), ylim=c(0,300), axes=F, ann=F)
axis(4, at=0:6*10, cex.axis=0.8)
mtext("Dose", 4, line=-1.5, at=30, cex=0.8)


