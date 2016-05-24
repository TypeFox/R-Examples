### R code from vignette source 'rgenoud.Rnw'

###################################################
### code chunk number 1: rgenoud.Rnw:488-495
###################################################
    claw <- function(xx) {
     x <- xx[1]
      y <- (0.46*(dnorm(x,-1.0,2.0/3.0) + dnorm(x,1.0,2.0/3.0)) +
      (1.0/300.0)*(dnorm(x,-0.5,.01) + dnorm(x,-1.0,.01) + dnorm(x,-1.5,.01)) +
     (7.0/300.0)*(dnorm(x,0.5,.07) + dnorm(x,1.0,.07) + dnorm(x,1.5,.07))) 
      return(y)
    }


###################################################
### code chunk number 2: rgenoud.Rnw:498-500
###################################################
library("rgenoud")
claw1   <- genoud(claw, nvars=1, max=TRUE, pop.size=3000)


###################################################
### code chunk number 3: rgenoud.Rnw:855-877 (eval = FALSE)
###################################################
## source("supplemental.R")
## gradcheck <- rep(TRUE,23)
## gradcheck[6:7] <- FALSE
## sizeset <- c(5000,10000)
## genset <- c(30,100)
## nreps <- 50
## gsarray <- array(NA, dim=c(length(sizeset), length(genset), 23, nreps))
## asc <- function(x) { as.character(x) }
## dimnames(gsarray) <- list(asc(sizeset), asc(genset), NULL, NULL);
## for (gsize in sizeset) {
##   for (ngens in genset) {
##     for (i in 1:23) {
##       for (j in 1:nreps) {
##         gsarray[as.character(gsize), as.character(ngens),i,j] <- 
##           genoud(testfuncs[[i]], nvars=testNparms[i], pop.size=gsize,
##                max.gen=ngens, hard.gen=TRUE, Domains=testbounds[[i]],
##                solution.tol=1e-6, boundary=1, gradient.check=gradcheck[i],
##               print=0)$value
##       }
##     }
##   }
## }


###################################################
### code chunk number 4: rgenoud.Rnw:1020-1027 (eval = FALSE)
###################################################
## LQDxmpl <- function(b) {
##   logistic <- function(x) { 1/(1+exp(-x)) }
##   sIQR <- function(y, yhat, n) {
##     IQR((y-yhat)/sqrt(yhat*(n-yhat)), na.rm=TRUE)
##   }
##  sIQR(y, m*logistic(x %*% b), m)
## }


###################################################
### code chunk number 5: rgenoud.Rnw:1032-1038 (eval = FALSE)
###################################################
## m <- 100
## x <- cbind(1,rnorm(1000),rnorm(1000))
## b1 <- c(.5, 1, -1)
## b2 <- c(0, -1, 1)
## logistic <- function(x) { 1/(1+exp(-x)) }
## y <- rbinom(1000, m, logistic(c(x[1:900,] %*% b1, x[901:1000,] %*% b2)))


###################################################
### code chunk number 6: rgenoud.Rnw:1047-1048 (eval = FALSE)
###################################################
## summary(glm1 <- glm(cbind(y,m-y) ~ x[,2] + x[,3], family="binomial"))


###################################################
### code chunk number 7: rgenoud.Rnw:1070-1073 (eval = FALSE)
###################################################
## suby <- y[1:900]
## subx <- x[1:900,]
## summary(glm2 <- glm(cbind(suby,m-suby) ~ subx[,2] + subx[,3], family="binomial"))


###################################################
### code chunk number 8: rgenoud.Rnw:1115-1132 (eval = FALSE)
###################################################
## dLQDxmpl <- function(b) {
##  eps <- 1e-10
##   logistic <- function(x) { 1/(1+exp(-x)) }
##   sIQR <- function(y, yhat, n) {
##     IQR((y-yhat)/sqrt(yhat*(n-yhat)), na.rm=TRUE)
##   }
##  dsIQR <- vector()
##   for (i in 1:length(b)) {
##    beps <- b
##     beps[i] <- b[i]+eps
##     dsIQR <-
##       c(dsIQR,
##         (sIQR(y, m*logistic(x %*% beps), m)-
##           sIQR(y, m*logistic(x %*% b), m))/eps)
##   }
##   return(dsIQR)
## }


###################################################
### code chunk number 9: rgenoud.Rnw:1136-1151 (eval = FALSE)
###################################################
## blen <- 3
## lenbseq <- length(bseq <- seq(-2,2,length=200))
## bseq3 <- seq(-1.2,-.9,length=200)
## bseq2 <- seq(.89,1.1,length=200)
## IQRarr <- IQRarrA <- array(NA, dim=c((1+blen), lenbseq, lenbseq))
## dimnames(IQRarrA) <- list(NULL, as.character(bseq), as.character(bseq))
## dimnames(IQRarr) <- list(NULL, as.character(bseq2), as.character(bseq3))
## for (i in 1:lenbseq) {
##   for (j in 1:lenbseq) {
##    IQRarrA[1,i,j] <- LQDxmpl(c(.5, bseq[i], bseq[j]))
##     IQRarrA[-1,i,j] <- dLQDxmpl(c(.5, bseq[i], bseq[j]))
##     IQRarr[1,i,j] <- LQDxmpl(c(.5, bseq2[i], bseq3[j]))
##     IQRarr[-1,i,j] <- dLQDxmpl(c(.5, bseq2[i], bseq3[j]))
##   }
## }


###################################################
### code chunk number 10: rgenoud.Rnw:1154-1165 (eval = FALSE)
###################################################
## par(mfrow=c(2,2), lwd=.1)
## contour(bseq,bseq, IQRarrA[1,,], main="IQR", xlab="b[2]", ylab="b[3]")
## contour(bseq,bseq, IQRarrA[3,,], main="partial derivative w/r/t b[2]",
##         xlab="b[2]", ylab="b[3]")
## loc2 <- (150:160)-5
## loc3 <- (40:50)+5
## contour(bseq[loc2],bseq[loc3], IQRarrA[3,loc2,loc3],
##         main="partial derivative w/r/t b[2]",
##         xlab="b[2]", ylab="b[3]")
## contour(bseq2,bseq3, IQRarr[3,,], main="partial derivative w/r/t b[2]",
##         xlab="b[2]", ylab="b[3]")


###################################################
### code chunk number 11: rgenoud.Rnw:1221-1224 (eval = FALSE)
###################################################
## LQD1  <-
##   genoud(LQDxmpl, nvars=3, max=FALSE, pop.size=2000, max.generations=300,
##         wait.generations=100, gradient.check=FALSE, print=1)


###################################################
### code chunk number 12: rgenoud.Rnw:1292-1295 (eval = FALSE)
###################################################
## LQD1  <-
##   genoud(LQDxmpl, nvars=3, max=FALSE, pop.size=10000, max.generations=1000,
##         wait.generations=300, gradient.check=FALSE, print=1)


