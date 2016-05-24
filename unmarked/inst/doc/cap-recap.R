### R code from vignette source 'cap-recap.Rnw'

###################################################
### code chunk number 1: cap-recap.Rnw:1-4
###################################################
options(width=70)
options(continue=" ")
library(tools)


###################################################
### code chunk number 2: cap-recap.Rnw:213-215
###################################################
alfl <- read.csv(system.file("csv", "alfl.csv", package="unmarked"))
head(alfl, 5)


###################################################
### code chunk number 3: cap-recap.Rnw:226-229
###################################################
alfl.covs <- read.csv(system.file("csv", "alflCovs.csv",
    package="unmarked"), row.names=1)
head(alfl.covs)


###################################################
### code chunk number 4: cap-recap.Rnw:243-249
###################################################
alfl$captureHistory <- paste(alfl$interval1, alfl$interval2, alfl$interval3, sep="")
alfl$captureHistory <- factor(alfl$captureHistory,
    levels=c("001", "010", "011", "100", "101", "110", "111"))
## Don't do this:
#levels(alfl$id) <- rownames(alfl.covs)
alfl$id <- factor(alfl$id, levels=rownames(alfl.covs))


###################################################
### code chunk number 5: cap-recap.Rnw:261-264
###################################################
alfl.v1 <- alfl[alfl$survey==1,]
alfl.H1 <- table(alfl.v1$id, alfl.v1$captureHistory)
head(alfl.H1, 5)


###################################################
### code chunk number 6: cap-recap.Rnw:298-310
###################################################
crPiFun <- function(p) {
    p1 <- p[,1]
    p2 <- p[,2]
    p3 <- p[,3]
    cbind("001" = (1-p1) * (1-p2) * p3,
          "010" = (1-p1) * p2     * (1-p3),
          "011" = (1-p1) * p2     * p3,
          "100" = p1     * (1-p2) * (1-p3),
          "101" = p1     * (1-p2) * p3,
          "110" = p1     * p2     * (1-p3),
          "111" = p1     * p2     * p3)
}


###################################################
### code chunk number 7: cap-recap.Rnw:319-322
###################################################
p <- matrix(0.4, 2, 3)
crPiFun(p)
rowSums(crPiFun(p))


###################################################
### code chunk number 8: cap-recap.Rnw:341-342
###################################################
o2y <- matrix(1, 3, 7)


###################################################
### code chunk number 9: cap-recap.Rnw:351-358
###################################################
library(unmarked)
intervalMat <- matrix(c('1','2','3'), 50, 3, byrow=TRUE)
class(alfl.H1) <- "matrix"
umf.cr1 <- unmarkedFrameMPois(y=alfl.H1,
    siteCovs=alfl.covs[,c("woody", "struct", "time.1", "date.1")],
    obsCovs=list(interval=intervalMat),
    obsToY=o2y, piFun="crPiFun")


###################################################
### code chunk number 10: cap-recap.Rnw:375-378
###################################################
M0 <- multinomPois(~1 ~1, umf.cr1)
Mt <- multinomPois(~interval-1 ~1, umf.cr1)
Mx <- multinomPois(~time.1 ~1, umf.cr1)


###################################################
### code chunk number 11: cap-recap.Rnw:387-388
###################################################
(M0.woody <- multinomPois(~1 ~woody, umf.cr1))


###################################################
### code chunk number 12: woody
###################################################
nd <- data.frame(woody=seq(0, 0.8, length=50))
E.abundance <- predict(M0.woody, type="state", newdata=nd, appendData=TRUE)
plot(Predicted ~ woody, E.abundance, type="l", ylim=c(0, 6),
     ylab="Alder flycatchers / plot", xlab="Woody vegetation cover")
lines(lower ~ woody, E.abundance, col=gray(0.7))
lines(upper ~ woody, E.abundance, col=gray(0.7))


###################################################
### code chunk number 13: cap-recap.Rnw:417-418
###################################################
backTransform(M0.woody, type="det")


###################################################
### code chunk number 14: cap-recap.Rnw:426-427
###################################################
round(getP(M0.woody), 2)[1,]


###################################################
### code chunk number 15: cap-recap.Rnw:454-465
###################################################
crPiFun.Mb <- function(p) { # p should have 3 columns
    pNaive <- p[,1]
    pWise <- p[,3]
    cbind("001" = (1-pNaive) * (1-pNaive) * pNaive,
          "010" = (1-pNaive) * pNaive     * (1-pWise),
          "011" = (1-pNaive) * pNaive     * pWise,
          "100" = pNaive     * (1-pWise)  * (1-pWise),
          "101" = pNaive     * (1-pWise)  * pWise,
          "110" = pNaive     * pWise      * (1-pWise),
          "111" = pNaive     * pWise      * pWise)
}


###################################################
### code chunk number 16: cap-recap.Rnw:478-484
###################################################
behavior <- matrix(c('Naive','Naive','Wise'), 50, 3, byrow=TRUE)
umf.cr1Mb <- unmarkedFrameMPois(y=alfl.H1,
    siteCovs=alfl.covs[,c("woody", "struct", "time.1")],
    obsCovs=list(behavior=behavior),
    obsToY=o2y, piFun="crPiFun.Mb")
M0 <- multinomPois(~1 ~1, umf.cr1Mb)


###################################################
### code chunk number 17: cap-recap.Rnw:489-490
###################################################
(Mb <- multinomPois(~behavior-1 ~1, umf.cr1Mb))


###################################################
### code chunk number 18: cap-recap.Rnw:496-497
###################################################
plogis(confint(Mb, type="det", method="profile"))


###################################################
### code chunk number 19: cap-recap.Rnw:544-576
###################################################
MhPiFun <- function(p) {
mu <- qlogis(p[,1]) # logit(p)
sig <- exp(qlogis(p[1,2]))
J <- ncol(p)
M <- nrow(p)
il <- matrix(NA, nrow=M, ncol=7)
dimnames(il) <- list(NULL, c("001","010","011","100","101","110","111"))
for(i in 1:M) {
  il[i,1] <- integrate( function(x) {
    (1-plogis(mu[i]+x))*(1-plogis(mu[i]+x))*plogis(mu[i]+x)*dnorm(x,0,sig)
  }, lower=-Inf, upper=Inf, stop.on.error=FALSE)$value
  il[i,2] <- integrate( function(x) {
    (1-plogis(mu[i]+x))*plogis(mu[i]+x)*(1-plogis(mu[i]+x))*dnorm(x,0,sig)
  }, lower=-Inf, upper=Inf, stop.on.error=FALSE)$value
  il[i,3] <- integrate( function(x) {
    (1-plogis(mu[i]+x))*plogis(mu[i]+x)*plogis(mu[i]+x)*dnorm(x,0,sig)
  }, lower=-Inf, upper=Inf, stop.on.error=FALSE)$value
  il[i,4] <- integrate( function(x) {
    plogis(mu[i]+x)*(1-plogis(mu[i]+x))*(1-plogis(mu[i]+x))*dnorm(x,0,sig)
  }, lower=-Inf, upper=Inf, stop.on.error=FALSE)$value
  il[i,5] <- integrate( function(x) {
    plogis(mu[i]+x)*(1-plogis(mu[i]+x))*plogis(mu[i]+x)*dnorm(x,0,sig)
  }, lower=-Inf, upper=Inf, stop.on.error=FALSE)$value
  il[i,6] <- integrate( function(x) {
    plogis(mu[i]+x)*plogis(mu[i]+x)*(1-plogis(mu[i]+x))*dnorm(x,0,sig)
  }, lower=-Inf, upper=Inf, stop.on.error=FALSE)$value
  il[i,7] <- integrate( function(x) {
    plogis(mu[i]+x)*plogis(mu[i]+x)*plogis(mu[i]+x)*dnorm(x,0,sig)
  }, lower=-Inf, upper=Inf, stop.on.error=FALSE)$value
}
return(il)
}


###################################################
### code chunk number 20: cap-recap.Rnw:694-703
###################################################
alfl.H <- table(alfl$id, alfl$captureHistory, alfl$survey)
alfl.Hmat <- cbind(alfl.H[,,1], alfl.H[,,2], alfl.H[,,3])
nVisits <- 3
o2yGMM <- kronecker(diag(nVisits), o2y)
umf.cr <- unmarkedFrameGMM(y=alfl.Hmat,
    siteCovs=alfl.covs[,c("woody", "struct")],
    yearlySiteCovs=list(date=alfl.covs[,3:5], time=alfl.covs[,6:8]),
    obsCovs=list(interval=cbind(intervalMat,intervalMat,intervalMat)),
    obsToY=o2yGMM, piFun="crPiFun", numPrimary=nVisits)


###################################################
### code chunk number 21: cap-recap.Rnw:722-723
###################################################
(fm1 <- gmultmix(~woody, ~1, ~time+date, umf.cr))


