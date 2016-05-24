### R code from vignette source 'nsRFA_ex01.Rnw'

###################################################
### code chunk number 1: nsRFA_ex01.Rnw:70-71
###################################################
library(nsRFA)


###################################################
### code chunk number 2: nsRFA_ex01.Rnw:81-82
###################################################
data(hydroSIMN)


###################################################
### code chunk number 3: nsRFA_ex01.Rnw:85-87 (eval = FALSE)
###################################################
## ls()
## help(hydroSIMN)


###################################################
### code chunk number 4: nsRFA_ex01.Rnw:107-111
###################################################
Dm <- parameters[,"Dm"]
logDm <- log(Dm)
sqrtDm <- sqrt(Dm)
sqrt3Dm <- Dm^(1/3)


###################################################
### code chunk number 5: nsRFA_ex01.Rnw:114-118
###################################################
attributes <- parameters[,-c(1,2)]
logattributes <- log(attributes[,-c(7:9)])
mixedattributes <- cbind(attributes, logattributes[,1])
 names(mixedattributes) <- c(names(attributes), "lnAm")


###################################################
### code chunk number 6: nontrasfregr
###################################################
nontrasfregr <- bestlm(Dm, mixedattributes, kmax=3, nbest=4); nontrasfregr


###################################################
### code chunk number 7: nsRFA_ex01.Rnw:127-136
###################################################
nregr <- dim(nontrasfregr$subselect)[1]
diagn <- data.frame(matrix(NA, nrow=nregr, ncol=2)); names(diagn) <- c("RMSE","RMSEjk")
for (i in 1:nregr){
 f <- paste("Dm ~", paste(colnames(nontrasfregr$subselect)[nontrasfregr$subselect[i,]], collapse=" + "))
 regr <- lm(f, mixedattributes)
 diagn[i,1] <- RMSE.lm(regr)
 diagn[i,2] <- RMSEjk.lm(regr)
}
diagn


###################################################
### code chunk number 8: multregr
###################################################
multregr <- bestlm(logDm, logattributes, kmax=3, nbest=4); multregr


###################################################
### code chunk number 9: nsRFA_ex01.Rnw:144-156
###################################################
nregr <- dim(multregr$subselect)[1]
diagn <- data.frame(matrix(NA, nrow=nregr, ncol=2)); names(diagn) <- c("RMSE","RMSEjk")
for (i in 1:nregr){
 f <- paste("logDm ~", paste(colnames(multregr$subselect)[multregr$subselect[i,]], 
                             collapse=" + "))
 regr <- lm(f, logattributes)
 fitt <- regr$fitted.values
 crossval <- jackknife1.lm(regr)
 diagn[i,1] <- RMSE(Dm, exp(fitt))
 diagn[i,2] <- RMSE(Dm, exp(crossval))
}
diagn


###################################################
### code chunk number 10: trasfregr_log
###################################################
trasfregr_log <- bestlm(logDm, mixedattributes, kmax=3, nbest=4); trasfregr_log


###################################################
### code chunk number 11: nsRFA_ex01.Rnw:163-175
###################################################
nregr <- dim(trasfregr_log$subselect)[1]
diagn <- data.frame(matrix(NA, nrow=nregr, ncol=2)); names(diagn) <- c("RMSE","RMSEjk")
for (i in 1:nregr){
 f <- paste("logDm ~", paste(colnames(trasfregr_log$subselect)[trasfregr_log$subselect[i,]], 
                             collapse=" + "))
 regr <- lm(f, mixedattributes)
 fitt <- regr$fitted.values
 crossval <- jackknife1.lm(regr)
 diagn[i,1] <- RMSE(Dm, exp(fitt))
 diagn[i,2] <- RMSE(Dm, exp(crossval))
}
diagn


###################################################
### code chunk number 12: trasfregr_sqrt
###################################################
trasfregr_sqrt <- bestlm(sqrtDm, mixedattributes, kmax=3, nbest=4); trasfregr_sqrt


###################################################
### code chunk number 13: nsRFA_ex01.Rnw:181-193
###################################################
nregr <- dim(trasfregr_sqrt$subselect)[1]
diagn <- data.frame(matrix(NA, nrow=nregr, ncol=2)); names(diagn) <- c("RMSE","RMSEjk")
for (i in 1:nregr){
 f <- paste("sqrtDm ~", paste(colnames(trasfregr_sqrt$subselect)[trasfregr_sqrt$subselect[i,]], 
                              collapse=" + "))
 regr <- lm(f, mixedattributes)
 fitt <- regr$fitted.values
 crossval <- jackknife1.lm(regr)
 diagn[i,1] <- RMSE(Dm, fitt^2)
 diagn[i,2] <- RMSE(Dm, crossval^2)
}
diagn


###################################################
### code chunk number 14: trasfregr_sqrt3
###################################################
trasfregr_sqrt3 <- bestlm(sqrt3Dm, mixedattributes, kmax=3, nbest=4); trasfregr_sqrt3


###################################################
### code chunk number 15: nsRFA_ex01.Rnw:199-211
###################################################
nregr <- dim(trasfregr_sqrt3$subselect)[1]
diagn <- data.frame(matrix(NA, nrow=nregr, ncol=2)); names(diagn) <- c("RMSE","RMSEjk")
for (i in 1:nregr){
 f <- paste("sqrt3Dm ~", paste(colnames(trasfregr_sqrt3$subselect)[trasfregr_sqrt3$subselect[i,]], 
                               collapse=" + "))
 regr <- lm(f, mixedattributes)
 fitt <- regr$fitted.values
 crossval <- jackknife1.lm(regr)
 diagn[i,1] <- RMSE(Dm, fitt^3)
 diagn[i,2] <- RMSE(Dm, crossval^3)
}
diagn


###################################################
### code chunk number 16: nsRFA_ex01.Rnw:217-218
###################################################
bestregr <- lm(sqrt3Dm ~ S2000 + IT + lnAm, mixedattributes); bestregr


###################################################
### code chunk number 17: nsRFA_ex01.Rnw:220-221
###################################################
summary(bestregr)


###################################################
### code chunk number 18: nsRFA_ex01.Rnw:225-227
###################################################
vif.lm(bestregr)
cor(bestregr$model[-1])


###################################################
### code chunk number 19: nsRFA_ex01.Rnw:230-231
###################################################
prt.lm(bestregr)


###################################################
### code chunk number 20: bestregression
###################################################
bestregr <- lm(logDm ~ Hm + NORD + IB, mixedattributes)
bestregr
summary(bestregr)


###################################################
### code chunk number 21: nsRFA_ex01.Rnw:241-242
###################################################
prt.lm(bestregr)


###################################################
### code chunk number 22: nsRFA_ex01.Rnw:244-246
###################################################
vif.lm(bestregr)
cor(bestregr$model[-1])


###################################################
### code chunk number 23: nsRFA_ex01.Rnw:249-250
###################################################
p_norm <- A2_GOFlaio(bestregr$residuals, dist="NORM"); p_norm


###################################################
### code chunk number 24: nsRFA_ex01.Rnw:253-257
###################################################
rmse <- RMSE(Dm, exp(bestregr$fitted.values))

predicted <- jackknife1.lm(bestregr)
rmse_jk <- RMSE(Dm, exp(predicted))


###################################################
### code chunk number 25: nsRFA_ex01.Rnw:264-280
###################################################
op <- par(mfrow=c(2,2))
 plot(bestregr$fitted.values, bestregr$residuals, xlab="Fitted", ylab="Residuals")
  abline(0,0,lty=3)

 normplot(bestregr$residuals, xlab="Residuals")
 
 plot(parameters[,c("Dm")], exp(bestregr$fitted.values), xlab="Originals", ylab="Fitted")
  abline(0,1,lty=3)
  intervals <- predinterval.lm(bestregr)
  intervals <- intervals[order(intervals[,1]),]
 
 plot(parameters[,c("Dm")], exp(predicted), xlab="Originals", ylab="Predicted")
  abline(0,1,lty=3)
  lines(exp(intervals[,c(1,2)]),lty=2)
  lines(exp(intervals[,c(1,3)]),lty=2)
par(op)


###################################################
### code chunk number 26: nsRFA_ex01.Rnw:303-306
###################################################
D <- annualflows["dato"][,]
y <- annualflows["anno"][,]
cod <- annualflows["cod"][,]


###################################################
### code chunk number 27: nsRFA_ex01.Rnw:313-314
###################################################
consistencyplot(y,cod)


###################################################
### code chunk number 28: nsRFA_ex01.Rnw:326-331
###################################################
ni <- tapply(D, cod, length)
annualflows15 <- annualflows[unsplit(ni, cod)>=15,]
parameters15 <- parameters[ni>=15,]
D15 <- annualflows15["dato"][,]
cod15 <- annualflows15["cod"][,]


###################################################
### code chunk number 29: nsRFA_ex01.Rnw:335-336
###################################################
LM15 <- data.frame(t(sapply(split(D15, cod15), Lmoments)))


###################################################
### code chunk number 30: nsRFA_ex01.Rnw:342-343
###################################################
plot(LM15[3:5])


###################################################
### code chunk number 31: nsRFA_ex01.Rnw:359-361
###################################################
Lspace.HWvsAD()
points(LM15[,4:3])


###################################################
### code chunk number 32: nsRFA_ex01.Rnw:374-375
###################################################
set.seed(10)


###################################################
### code chunk number 33: nsRFA_ex01.Rnw:377-379
###################################################
D15adim <- D15/unsplit(tapply(D15, cod15, mean), cod15)
HWs <- HW.tests(D15adim, cod15)[1]; HWs


###################################################
### code chunk number 34: nsRFA_ex01.Rnw:383-384
###################################################
bestlm(LM15[,"lcv"], parameters15[,3:16], kmax=3)


###################################################
### code chunk number 35: nsRFA_ex01.Rnw:387-389
###################################################
bestlm(as.numeric(AD.dist(D15,cod15)), data.frame(apply(parameters15[,3:16], 2, dist)), 
       kmax=3)


###################################################
### code chunk number 36: mantel.test
###################################################
Y <- AD.dist(D15,cod15)
X <- data.frame(apply(parameters15[,c("Hm","Ybar")],2,dist))
datamantel <- cbind(as.numeric(Y),X)
regrmantel <- lm(Y ~ Hm + Ybar, datamantel)
#summary(regrmantel)
mantel.lm(regrmantel, Nperm=100)


###################################################
### code chunk number 37: clusteranalysis
###################################################
param <- parameters15[c("Hm","Ybar")]
n <- dim(param)[1]; k <- dim(param)[2]
param.norm <- (param - matrix(sapply(param, mean), nrow=n, ncol=k, byrow=TRUE))/
              matrix(sapply(param, sd), nrow=n, ncol=k, byrow=TRUE)


###################################################
### code chunk number 38: nsRFA_ex01.Rnw:410-411
###################################################
set.seed(10)


###################################################
### code chunk number 39: nsRFA_ex01.Rnw:413-424
###################################################
nclusters=1
while (max(HWs) > 2.1) {
 nclusters <- nclusters+1
 clusters <- traceWminim(param.norm, nclusters)
 indclusters <- unsplit(clusters, cod15)
 HWs <- rep(NA, nclusters)
 for (i in unique(clusters)) {
  HWs[i] <- HW.tests(D15adim[indclusters==i], cod15[indclusters==i])[1]
 }
 print(HWs)
} 


###################################################
### code chunk number 40: nsRFA_ex01.Rnw:431-433
###################################################
regLM15 <- t(sapply(split(D15adim, indclusters), Lmoments))[,3:5]
regLM15


###################################################
### code chunk number 41: nsRFA_ex01.Rnw:437-440
###################################################
for (i in 1:nclusters) {
 print(regionalLmoments(D15adim[indclusters==i], cod15[indclusters==i])[3:5])
}


###################################################
### code chunk number 42: nsRFA_ex01.Rnw:448-472
###################################################
op <- par(mfrow=c(2,2))
 plot(parameters15[c("Hm","Ybar")], col=clusters, pch=clusters, cex=0.6,
      main="Clusters in the space of classification variables", cex.main=1, font.main=1)
  grid()
  points(tapply(parameters15["Hm"][,], clusters, mean), 
         tapply(parameters15["Ybar"][,], clusters, mean),
         col=c(1:nclusters), pch=c(1:nclusters))
 legend("topleft", paste("clust ",c(1:nclusters)), 
        col=c(1:nclusters), pch=c(1:nclusters), bty="n")

 plot(parameters15[c("Xbar","Ybar")], col=clusters, pch=clusters, cex=0.6,
      main="Clusters in geographical space", cex.main=1, font.main=1)
  grid()

 plot(LM15[,4:3], pch=clusters, col=clusters, cex=0.6,
      main="Clusters in L-moments space", cex.main=1, font.main=1)
  points(regLM15[,2:1], col=c(1:nclusters), pch=c(1:nclusters))
  grid()

 plot(LM15[,4:5], pch=clusters, col=clusters, cex=0.6,
      main="Clusters in L-moments space", cex.main=1, font.main=1)
  points(regLM15[,2:3], col=c(1:nclusters), pch=c(1:nclusters))
  grid()
par(op)


###################################################
### code chunk number 43: nsRFA_ex01.Rnw:489-493
###################################################
Lmoment.ratio.diagram()
 points(regLM15[,2:3], col=c(1:nclusters), pch=c(1:nclusters))
 legend("bottomleft",paste("clust ", c(1:nclusters)), 
        col=c(1:nclusters), pch=c(1:nclusters), bty="n")


###################################################
### code chunk number 44: GOFlaio2004
###################################################
for (i in 1:nclusters) {
  GOFA2_P3 <- A2_GOFlaio(D15adim[indclusters==i], dist="P3")
  cat(paste("\np(A2) for Cluster ", i, ":\n", sep=""))
  print(GOFA2_P3)
}


###################################################
### code chunk number 45: nsRFA_ex01.Rnw:538-544
###################################################
paramgamma=NULL
for (i in 1:nclusters) {
 paramgamma[[i]] <- par.gamma(1, regLM15[i,1], regLM15[i,2])
 cat(paste("\nCluster",i,":\n"))
 print(format(paramgamma[[i]][1:3]))
}


###################################################
### code chunk number 46: nsRFA_ex01.Rnw:547-552
###################################################
for (i in 1:nclusters) {
 cat(paste("\nCluster",i,":\n"))
 print(format(par2mom.gamma(paramgamma[[i]]$xi, 
       paramgamma[[i]]$beta, paramgamma[[i]]$alfa)))
}


###################################################
### code chunk number 47: nsRFA_ex01.Rnw:560-573
###################################################
op <- par(mfrow=c(2,2))
 for (i in 1:nclusters) {
  FF <- F.gamma(D15adim[indclusters==i], paramgamma[[i]]$xi, 
                paramgamma[[i]]$beta, paramgamma[[i]]$alfa)
  regionalplotpos(D15adim[indclusters==i], cod15[indclusters==i], 
                  xlab=paste("cluster", i),
                  main="Empirical distributions", cex.main=1, font.main=1)
  lines(sort(D15adim[indclusters==i]), sort(FF))
  nomi <- names(clusters)[clusters==i]
  legend("bottomright", legend=nomi, pch=c(1:length(nomi)), 
         col=c(1:length(nomi)), bty="n", cex=.9)
 }
par(op)


###################################################
### code chunk number 48: nsRFA_ex01.Rnw:588-601
###################################################
op <- par(mfrow=c(2,2))
 for (i in 1:nclusters) {
  Fs <- seq(0.001,0.999,by=.001)
  regionalnormplot(D15adim[indclusters==i], cod15[indclusters==i], 
                   xlab=paste("cluster", i),
                   main="Empirical distributions", cex.main=1, font.main=1)
  normpoints(invF.gamma(Fs, paramgamma[[i]]$xi, paramgamma[[i]]$beta, 
                        paramgamma[[i]]$alfa), type="l")
  nomi <- names(clusters)[clusters==i]
  legend("bottomright", legend=nomi, pch=c(1:length(nomi)), 
         col=c(1:length(nomi)), bty="n", cex=.9)
 }
par(op)


###################################################
### code chunk number 49: nsRFA_ex01.Rnw:616-626
###################################################
spess=c(1, 1.5, 2, 1.3)
Fs <- seq(0.001,0.999,by=.001)
lognormplot(D15adim, line=FALSE, type="n", )
for (i in 1:nclusters) {
 qq <- invF.gamma(Fs, paramgamma[[i]]$xi, paramgamma[[i]]$beta, 
                  paramgamma[[i]]$alfa)
 normpoints(qq, type="l", lty=i, col=i, lwd=spess[i])
}
legend("bottomright", paste("cluster ", c(1:nclusters)), 
       col=c(1:nclusters), lty=c(1:nclusters), lwd=spess, bty="n")


