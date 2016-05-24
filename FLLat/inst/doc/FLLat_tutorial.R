### R code from vignette source 'FLLat_tutorial.rnw'

###################################################
### code chunk number 1: FLLat_tutorial.rnw:110-114
###################################################
library(FLLat)
data(simaCGH)
result.pve <- FLLat.PVE(simaCGH,J.seq=1:ncol(simaCGH))
plot(result.pve)


###################################################
### code chunk number 2: FLLat_tutorial.rnw:139-142
###################################################
result.bic <- FLLat.BIC(simaCGH,J=5)
result.bic$lam1
result.bic$lam2


###################################################
### code chunk number 3: FLLat_tutorial.rnw:153-154
###################################################
plot(result.bic$opt.FLLat)


###################################################
### code chunk number 4: FLLat_tutorial.rnw:164-165
###################################################
plot(result.bic$opt.FLLat,type="weights")


###################################################
### code chunk number 5: FLLat_tutorial.rnw:227-230
###################################################
result.fdr <- FLLat.FDR(simaCGH,result.bic$opt.FLLat)
result.fdr$thresh.control
plot(result.fdr)


###################################################
### code chunk number 6: FLLat_tutorial.rnw:267-273
###################################################
tr.dat <- simaCGH[,1:15]
tst.dat <- simaCGH[,16:20]
result.tr <- FLLat(tr.dat,J=5,lam1=1,lam2=9)
tst.pred <- predict(result.tr,newY=tst.dat)
plot(tst.dat[,1],xlab="Probe",ylab="Y")
lines(tst.pred$pred.Y[,1],col="red",lwd=3)


