### R code from vignette source 'laplace.Rnw'

###################################################
### code chunk number 1: laplace.Rnw:33-35
###################################################
options(continue="  ", width=60)
options(SweaveHooks=list(fig=function() par(mar=c(5.1, 4.1, .3, 1.1))))


###################################################
### code chunk number 2: laplace.Rnw:141-167
###################################################
library(coxme)
set.seed(1953)  # an auspicious birth year :-)
mkdata <- function(n, beta=c(.4, .1), sitehaz=c(.5,1.5, 2,1)) {
    nsite <- length(sitehaz)
    site <- rep(1:nsite, each=n)
    trt1 <- rep(0:1, length=n*nsite)
    hazard <- sitehaz[site] + beta[1]*trt1 + beta[2]*trt1 * (site-mean(site))
    stime <- rexp(n*nsite, exp(hazard))
    q80 <- quantile(stime, .8)
    data.frame(site=site,
               trt = trt1,
               futime= pmin(stime, q80),
               status= ifelse(stime>q80, 0, 1),
               hazard=hazard
               )
}
trdata <- mkdata(150)  #150 enrolled per site
fit1 <- coxme(Surv(futime, status) ~ trt + (1| site/trt), trdata)
print(fit1)

# Show the true and estimated per-site intercepts
true <- c(.5, 1.5, 2, 1) - mean(c(.5, 1.5, 2, 1))
bcoef <- ranef(fit1)[[2]]
temp <- rbind(true, bcoef)
dimnames(temp) <- list(c("True", "Estimated"), paste("Site",1:4))
round(temp,2)


###################################################
### code chunk number 3: fig1
###################################################
getOption("SweaveHooks")[["fig"]]()
xx <- seq(-1, 1, length=101) #vary b from -1 to 1
profile <- matrix(0, nrow=101, ncol=8) #to store curves
bcoef <- unlist(ranef(fit1))
indx <- -1 + trdata$trt + 2*trdata$site  #random treatment effect index
Ainv <- diag(1/rep(unlist(VarCorr(fit1)), c(8,4)))
for (i in 1:4) {
    tcoef <- bcoef
    for (j in 1:101) {
        tcoef[i+8] <- xx[j]  #reset single coef
        eta <- fixef(fit1)*trdata$trt + tcoef[trdata$site+8] +
            tcoef[indx]
        tfit <- coxph(Surv(futime, status) ~ offset(eta), data= trdata)
        
        profile[j,i] <- tfit$loglik - .5*tcoef%*% Ainv %*% tcoef
        profile[j, i+4] <- fit1$loglik[3] - 
            .5*sum(((tcoef-bcoef) %*% fit1$hmat[1:12, 1:12])^2)
    }
}
matplot(xx, profile-fit1$loglik[3], type='l', lty=c(1,1,1,1,2,2,2,2), col=1:4,
        ylim=c(-40,0),
        xlab="b", ylab="LPPL - max")


###################################################
### code chunk number 4: laplace.Rnw:286-289
###################################################
fit1b <- coxme(Surv(futime, status) ~ trt + (1 | site/trt),
               data=trdata, refine.n=500)
fit1b$refine


###################################################
### code chunk number 5: laplace.Rnw:306-315
###################################################
efit2 <- coxme(Surv(y, uncens) ~ trt + (1|center), eortc,
                refine.n=100)
efit2$refine

efit3 <- coxme(Surv(y, uncens) ~ trt + (1|center/trt), eortc,
               refine.n=100)
efit3$refine

efit3


###################################################
### code chunk number 6: laplace.Rnw:342-346
###################################################
lfit1 <- coxph(Surv(time, status) ~ age + ph.ecog + wt.loss, lung)
lfit2 <- coxme(Surv(time, status) ~ age + (ph.ecog |1) +
               (wt.loss |1), data=lung, refine.n=100)
lfit2$refine


###################################################
### code chunk number 7: laplace.Rnw:352-356
###################################################
print(lfit2, rcoef=TRUE)

signif(rbind(coef(lfit1), 
             c(fixef(lfit2), unlist(ranef(lfit2)))),2)


###################################################
### code chunk number 8: cgd
###################################################
getOption("SweaveHooks")[["fig"]]()
cfit <- coxme(Surv(tstart, tstop, status) ~ treat + age +
                  (1 | id), data=cgd, refine.n=500, refine.detail=TRUE)
cfit$refine
2*(diff(cfit$loglik[1:2]))

temp <- cfit$refine.detail
e1 <- (temp$loglik - temp$penalty1) - cfit$loglik[2]
e2 <- (cfit$loglik[3] - temp$penalty2) - cfit$loglik[2]
ssqrt <- function(x) sign(x)*sqrt(abs(x))  #signed square root
plot(ssqrt(e1), ssqrt(e2), xlab="sqrt(e1)", ylab="sqrt(e2)")
abline(0,1)


###################################################
### code chunk number 9: cgd2
###################################################
getOption("SweaveHooks")[["fig"]]()
ss <- seq(.3, 1.3, length=25)
tmat <- matrix(0, nrow=25, ncol=3)
for (i in 1:25) {
    tfit <- coxme(Surv(tstart, tstop, status) ~ treat + age + (1|id),
                  cgd, vfixed=ss[i]^2, refine.n=1000)
    tmat[i,] <- c(diff(tfit$loglik[1:2]), tfit$refine)
}
temp1 <- tmat[,1] + tmat[,2]  #corrected IPL
temp2 <- tmat[,1] + tmat[,2] + cbind(-2*tmat[,3], 2*tmat[,3]) # .955 CI
matplot(ss, cbind(tmat[,1], temp1), pch='o', col=1:2,
        ylim=range(tmat[,1], temp2), 
        xlab="Std of random effect",
        ylab="IPL - Null")
segments(ss, temp2[,1], ss, temp2[,2], lty=2, col=2)
lines(smooth.spline(ss, temp1, df=5), col=2)
abline(h= diff(cfit$loglik[1:2]) - qchisq(.95, 1)/2, lty=2)


###################################################
### code chunk number 10: laplace.Rnw:470-473
###################################################
cfit1 <- coxph(Surv(time, status) ~ rx + nodes + extent +
         strata(etype) + cluster(id), colon)
cfit1


###################################################
### code chunk number 11: laplace.Rnw:486-494
###################################################
cfit2 <- coxme(Surv(time, status) ~ rx + nodes + extent + 
               strata(etype) + (1|id), colon,
               refine.n=500)
cfit2$refine

print(cfit2)

round(quantile(ranef(cfit2)[[1]], 0:8/8), 2)


