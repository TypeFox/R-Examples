# This is file ../spam/tests/jss_areal_counts.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     




# JSS article: 
#     "Pitfalls in the implementation of Bayesian
#      hierarchical modeling of areal count data.
#      An illustration using BYM and Leroux models."
#
# test the MCMC sampler from the paper with 30 iterations.


# SETUP:
library("spam")
spam.options(structurebased=TRUE)

# BYM ---------------------------------------------
data(Oral); attach(Oral) 
path <- system.file("demodata/germany.adjacency", package = "spam")
A <- adjacency.landkreis(path); n <- dim(A)[1]

set.seed(2)
hyperA <- c(1, 1); hyperB <- c(0.5, .01)  
totalg <- 30

upost <- vpost <- array(0, c(totalg, n))  
kpost <- array(NA, c(totalg, 2)); accept <- rep(NA, totalg)
upost[1,] <- vpost[1,] <- rep(.001, 544); kpost[1,] <- c(10, 100) 

eta <- upost[1,] + vpost[1,]
C <- exp(eta) * E; diagC <- diag.spam(c(rep(0, n), C))
b <- c( rep(0, n), Y + (eta - 1) * C)
Qu <- R <- precmat.IGMRFirreglat(A); pad(Qu) <- c(2 * n, 2 * n)
Qv <- as.spam(rbind(cbind( diag(n), -diag(n)),
                      cbind(-diag(n),  diag(n))))
Q <- kpost[1,1] * Qu + kpost[1,2] * Qv + diagC
struct <- chol(Q, memory = list(nnzcolindices = 6467))
uRuHalf <- t(upost[1,]) %*% (R %*% upost[1,]) / 2
vvHalf <- t(vpost[1,]) %*% vpost[1,] / 2
postshape <- hyperA + c(n - 1, n) / 2

for (i in 2:totalg) {
     kpost[i,] <- rgamma(2, postshape, hyperB + c(uRuHalf, vvHalf))
 
     etaTilde <- eta
     for(index in 1:2){
         C <- E * exp(etaTilde)
         diagC <- diag.spam(c(rep(0, n), C))
         b <- c(rep(0, 544), Y + (etaTilde - 1) * C)
         Q <- kpost[i,1] * Qu + kpost[i,2] * Qv + diagC
         etaTilde <- c(solve.spam(Q, b,
                                  Rstruct = struct))[1:n + n]
     }
 
     C <- exp(etaTilde) * E; diagC <- diag.spam(c(rep(0, n), C))
     b <- c( rep(0, n), Y + (etaTilde - 1) * C)
     Q <- kpost[i,1] * Qu + kpost[i,2] * Qv + diagC

     x_ <- c(rmvnorm.canonical(1, b, Q, Rstruct = struct))
     upost[i,] <- x_[1:n]; eta_ <- x_[1:n + n]; vpost[i,] <- eta_ - upost[i,]
     uRuHalf_ <- t(upost[i,]) %*% (R %*% upost[i,]) / 2
     vvHalf_ <- t(vpost[i,]) %*% vpost[i,] / 2

     etaTilde_ <- eta_
     for(index in 1:2){
       C_ <- E * exp(etaTilde_)
       diagC_ <- diag.spam(c(rep(0, n), C_))
       b_ <- c(rep(0, 544), Y + (etaTilde_ - 1) * C_)
       Q_<- kpost[i,1] * Qu + kpost[i,2] * Qv + diagC_
       etaTilde_ <- c(solve.spam(Q_, b_,
                                Rstruct = struct))[1:n + n]
     }
     
     C_ <- exp(etaTilde_) * E; diagC_ <- diag.spam(c(rep(0, n), C_))
     b_ <- c( rep(0, n), Y + (etaTilde_ - 1) * C_)
     Q_ <- kpost[i,1] * Qu + kpost[i,2] * Qv + diagC_

     logPost_ <- sum(Y * eta_ - E * exp(eta_)) -
         kpost[i,1] * uRuHalf_ - kpost[i, 2] * vvHalf_
     logPost  <- sum(Y * eta - E * exp(eta)) - kpost[i,1] * uRuHalf -
         kpost[i,2] * vvHalf
     logApproxX_ <- - kpost[i,1] * uRuHalf_ - kpost[i,2] * vvHalf_ -
         sum(.5 * eta_^2 * C) + sum(b * eta_)
     logApproxX  <- - kpost[i,1] * uRuHalf  - kpost[i,2] * vvHalf -
         sum(.5 * eta^2 * C_) + sum(b_ * eta)
     logAlpha <- min(0, logPost_ - logPost + logApproxX - logApproxX_)
 
     if (log(runif(1)) < logAlpha) {    
         uRuHalf <- uRuHalf_;  vvHalf <- vvHalf_
         eta <- eta_; b <- b_; C <- C_; accept[i] <- 1
     } else{                            
       accept[i] <- 0; upost[i,] <- upost[i-1,]; vpost[i,] <- vpost[i-1,]} 
 }
 
# values of  30th iteration
head(eta)
tail(b)
head(C)
tail(accept)
tail(upost[30,])
tail(vpost[30,])
sum(accept[-1])
sum(upost)

