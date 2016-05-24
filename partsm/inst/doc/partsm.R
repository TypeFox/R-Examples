### R code from vignette source 'partsm.Rnw'

###################################################
### code chunk number 1: partsm.Rnw:176-179
###################################################
library(partsm)
data("gergnp")
lgergnp <- log(gergnp, base=exp(1))


###################################################
### code chunk number 2: partsm.Rnw:190-200
###################################################
detcomp <- list(regular=c(0,0,0), seasonal=c(1,0), regvar=0)
aic <- bic <- Fnextp <- Fpval <- rep(NA, 4)
for(p in 1:4){
  lmpar <- fit.ar.par(wts=lgergnp, detcomp=detcomp, type="PAR", p=p)
  aic[p] <- AIC(lmpar@lm.par, k=2)
  bic[p] <- AIC(lmpar@lm.par, k=log(length(residuals(lmpar@lm.par))))
  Fout <- Fnextp.test(wts=lgergnp, detcomp=detcomp, p=p, type="PAR")
  Fnextp[p] <- Fout@Fstat
  Fpval[p] <- Fout@pval
}


###################################################
### code chunk number 3: partsm.Rnw:229-235
###################################################
dcsi <- list(regular=c(0,0,0), seasonal=c(1,0), regvar=0)
out.Fparsi  <- Fpar.test(wts=lgergnp, detcomp=dcsi, p=2)
show(out.Fparsi)
dcsit <- list(regular=c(0,0,0), seasonal=c(1,1), regvar=0)
out.Fparsit <- Fpar.test(wts=lgergnp, detcomp=dcsit, p=2)
show(out.Fparsit)


###################################################
### code chunk number 4: partsm.Rnw:244-247
###################################################
par2 <- fit.ar.par(wts=lgergnp, type="PAR", p=2, detcomp=detcomp)
Fsh.out <- Fsh.test(res=residuals(par2@lm.par), s=frequency(lgergnp))
show(Fsh.out)


###################################################
### code chunk number 5: partsm.Rnw:259-262
###################################################
out.par <- fit.ar.par(wts=lgergnp, type="PAR", detcomp=detcomp, p=2)
out.MV <- PAR.MVrepr(out.par)
out.MV


###################################################
### code chunk number 6: partsm.Rnw:284-286
###################################################
out.LR <- LRurpar.test(wts=lgergnp, detcomp=detcomp, p=2)
show(out.LR)


###################################################
### code chunk number 7: partsm.Rnw:291-293
###################################################
Fpari1.out <- Fpari.piar.test(wts=lgergnp, detcomp=detcomp, p=2, type="PARI1")
show(Fpari1.out)


###################################################
### code chunk number 8: partsm.Rnw:349-352
###################################################
out.piar <- fit.piar(wts=lgergnp, detcomp=detcomp, p=2)
out.MV <- PAR.MVrepr(out.piar)
out.MV


###################################################
### code chunk number 9: pdiffplot
###################################################
plotpdiff(out.piar)


###################################################
### code chunk number 10: partsm.Rnw:374-376
###################################################
out.pred <- predictpiar(wts=lgergnp, p=2, hpred=24)
show(out.pred)


###################################################
### code chunk number 11: partsm.Rnw:381-386
###################################################
out.pred@wts <- exp(1)^out.pred@wts
out.pred@fcast <- exp(1)^out.pred@fcast
out.pred@ucb <- exp(1)^out.pred@ucb
out.pred@lcb <- exp(1)^out.pred@lcb
plotpredpiar(out.pred)


###################################################
### code chunk number 12: predplot
###################################################
plotpredpiar(out.pred)


