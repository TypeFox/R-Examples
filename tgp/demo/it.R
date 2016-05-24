###################################################
### chunk number 1: 
###################################################
library(tgp)
library(maptree)
#options(width=65)
seed <- 0; set.seed(seed)


###################################################
### chunk number 2: 
###################################################
geo <- default.itemps(type="geometric")
har <- default.itemps(type="harmonic")
sig <- default.itemps(type="sigmoidal")


###################################################
### chunk number 3: it-itemps
###################################################
par(mfrow=c(2,1))
all <- cbind(geo$k, har$k, sig$k)
matplot(all, pch=21:23,
        main="inv-temp ladders", xlab="indx", ylab="itemp")
legend("topright", pch=21:23, 
       c("geometric","harmonic","sigmoidal"), col=1:3)
matplot(log(all), pch=21:23,
        main="log(inv-temp) ladders", xlab="indx", ylab="itemp")


###################################################
### chunk number 4: 
###################################################
graphics.off()


###################################################
### chunk number 5: 
###################################################
ESS <- function(w)
{
  mw <- mean(w)
  cv2 <- sum((w-mw)^2)/((length(w)-1)*mw^2)
  ess <- length(w)/(1+cv2)
  return(ess)
}


###################################################
### chunk number 6: 
###################################################
exp2d.data<-exp2d.rand() 
X<-exp2d.data$X 
Z<-exp2d.data$Z 


###################################################
### chunk number 7: 
###################################################
its <- default.itemps(m=10)
exp.btlm <- btlm(X=X,Z=Z, bprior="b0", R=2, itemps=its, pred.n=FALSE,
                 BTE=c(1000,3000,2)) 


###################################################
### chunk number 8: 
###################################################
exp.btlm$ess


###################################################
### chunk number 9: 
###################################################
library(MASS)
moto.it <- btgpllm(X=mcycle[,1], Z=mcycle[,2], BTE=c(2000,52000,10),
        bprior="b0", R=3, itemps=geo, trace=TRUE, pred.n=FALSE, verb=0)


###################################################
### chunk number 10: 
###################################################
moto.it$ess$combined


###################################################
### chunk number 11: 
###################################################
p <- moto.it$trace$post
ESS(p$wlambda)


###################################################
### chunk number 12: 
###################################################
ESS(p$w)


###################################################
### chunk number 13: 
###################################################
as.numeric(c(sum(p$itemp == 1), moto.it$ess$each[1,2:3]))


###################################################
### chunk number 14: 
###################################################
moto.reg <- btgpllm(X=mcycle[,1], Z=mcycle[,2], BTE=c(2000,52000,10),
        R=3, bprior="b0", trace=TRUE, pred.n=FALSE, verb=0)


###################################################
### chunk number 15: 
###################################################
L <- length(p$height)
hw <- suppressWarnings(sample(p$height, L, prob=p$wlambda, replace=TRUE))
b <- hist2bar(cbind(moto.reg$trace$post$height, p$height, hw))


###################################################
### chunk number 16: it-moto-height
###################################################
barplot(b, beside=TRUE, col=1:3, xlab="tree height", ylab="counts", 
         main="tree heights encountered")
legend("topright", c("reg MCMC", "All Temps", "IT"), fill=1:3)


###################################################
### chunk number 17: 
###################################################
graphics.off()


###################################################
### chunk number 18: it-moto-ktrace
###################################################
plot(log(moto.it$trace$post$itemp), type="l", ylab="log(k)", xlab="samples",
     main="trace of log(k)")


###################################################
### chunk number 19: 
###################################################
graphics.off()


###################################################
### chunk number 20: it-moto-khist
###################################################
b <- itemps.barplot(moto.it, plot.it=FALSE)
barplot(t(cbind(moto.it$itemps$counts, b)), col=1:2,
        beside=TRUE, ylab="counts", xlab="itemps", 
        main="inv-temp observation counts")
legend("topleft", c("observation counts", "posterior samples"), fill=1:2)


###################################################
### chunk number 21: 
###################################################
graphics.off()


###################################################
### chunk number 22: 
###################################################
moto.it.sig <- btgpllm(X=mcycle[,1], Z=mcycle[,2], BTE=c(2000,52000,10),
                      R=3, bprior="b0", krige=FALSE, itemps=sig, verb=0)


###################################################
### chunk number 23: 
###################################################
moto.it.sig$ess$combined


###################################################
### chunk number 24: it-moto-pred
###################################################
plot(moto.it.sig)


###################################################
### chunk number 25: 
###################################################
graphics.off()


###################################################
### chunk number 26: 
###################################################
Xcand <- lhs(10000, rbind(c(-6,6),c(-6,6)))
X <- dopt.gp(400, X=NULL, Xcand)$XX
Z <- exp2d.Z(X)$Z


###################################################
### chunk number 27: 
###################################################
exp.reg <- btgpllm(X=X, Z=Z, BTE=c(2000,52000,10), bprior="b0", 
                   trace=TRUE, krige=FALSE, R=10, verb=0)


###################################################
### chunk number 28: it-exp-pred
###################################################
plot(exp.reg)


###################################################
### chunk number 29: 
###################################################
graphics.off()


###################################################
### chunk number 30: 
###################################################
h <- exp.reg$post$height[which.max(exp.reg$posts$lpost)]
h


###################################################
### chunk number 31: it-exp-mapt
###################################################
tgp.trees(exp.reg, "map")


###################################################
### chunk number 32: 
###################################################
graphics.off()


###################################################
### chunk number 33: 
###################################################
its <- default.itemps(k.min=0.02)
exp.it <- btgpllm(X=X, Z=Z, BTE=c(2000,52000,10), bprior="b0", 
               trace=TRUE, krige=FALSE, itemps=its, R=10, verb=0)


###################################################
### chunk number 34: 
###################################################
exp.it$gpcs
exp.reg$gpcs


###################################################
### chunk number 35: 
###################################################
p <- exp.it$trace$post
data.frame(ST=sum(p$itemp == 1), nIT=ESS(p$w), oIT=exp.it$ess$combined)


###################################################
### chunk number 36: 
###################################################
L <- length(p$height)
hw <- suppressWarnings(sample(p$height, L, prob=p$wlambda, replace=TRUE))
b <- hist2bar(cbind(exp.reg$trace$post$height, p$height, hw))


###################################################
### chunk number 37: it-exp-height
###################################################
barplot(b, beside=TRUE, col=1:3, xlab="tree height", ylab="counts", 
         main="tree heights encountered")
legend("topright", c("reg MCMC", "All Temps", "IT"), fill=1:3)


###################################################
### chunk number 38: 
###################################################
graphics.off()


###################################################
### chunk number 39: it-exp-trace-height
###################################################
ylim <- range(p$height, exp.reg$trace$post$height)
plot(p$height, type="l", main="trace of tree heights", 
     xlab="t", ylab="height", ylim=ylim)
lines(exp.reg$trace$post$height, col=2)
legend("topright", c("tempered", "reg MCMC"), lty=c(1,1), col=1:2)


###################################################
### chunk number 40: 
###################################################
graphics.off()


###################################################
### chunk number 41: it-expit-pred
###################################################
plot(exp.it)


###################################################
### chunk number 42: it-expit-trees
###################################################
tgp.trees(exp.it, "map")


###################################################
### chunk number 43: 
###################################################
graphics.off()


