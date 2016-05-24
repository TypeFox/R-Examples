# This is file ../spam/tests/demo_jss15-BYM.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     

options( echo=FALSE)
library( spam, warn.conflict=FALSE)



# This is the MCMC sampler presented in Section 2.2 of the article:
#
# Florian Gerber, Reinhard Furrer (2015). Pitfalls in the Implementation
# of Bayesian Hierarchical Modeling of Areal Count Data: An Illustration
# Using BYM and Leroux Models. Journal of Statistical Software, 
# Code Snippets, 63(1), 1–32. URL http://www.jstatsoft.org/v63/c01/.
#
# Note: For illustration we set
# number of generated samples: 5'000
# number of burnin samples:      500
# thinning:                       10
# This takes 1-2 minutes of computation time.
#
# In the JSS article we used:
# number of generated samples: 300'000
# number of burnin samples:     15'000
# thinning:                         20

# invisible(readline(prompt = "Type  <Return>\t to continue : "))

# SETUP:
spam.options(structurebased=TRUE)

# Besag-York-Molie model (BYM) 
# load data
data(Oral); attach(Oral) 
path <- system.file("demodata/germany.adjacency", package = "spam")
A <- adjacency.landkreis(path); n <- dim(A)[1]

set.seed(2)
# hyper priors
hyperA <- c(1, 1); hyperB <- c(0.5, .01)  

# sampler length, burnin and thinning
totalg <- 30
burnin <- 500
thin <- 10

# variable to store samples
upost <- vpost <- array(0, c(totalg, n))  
kpost <- array(NA, c(totalg, 2)); accept <- rep(NA, totalg)

# initial values
upost[1,] <- vpost[1,] <- rep(.001, 544); kpost[1,] <- c(10, 100) 

# precalculate some quantities
eta <- upost[1,] + vpost[1,]
C <- exp(eta) * E; diagC <- diag.spam(c(rep(0, n), C))
b <- c( rep(0, n), Y + (eta - 1) * C)
Qu <- R <- precmat.IGMRFirreglat(A); pad(Qu) <- c(2 * n, 2 * n)
Qv <- as.spam(rbind(cbind( diag(n), -diag(n)),
                      cbind(-diag(n),  diag(n))))
Q <- kpost[1,1] * Qu + kpost[1,2] * Qv + diagC

# store symbolic cholesky factorization of Q
struct <- chol(Q, memory = list(nnzcolindices = 6467))

uRuHalf <- t(upost[1,]) %*% (R %*% upost[1,]) / 2
vvHalf <- t(vpost[1,]) %*% vpost[1,] / 2
postshape <- hyperA + c(n - 1, n) / 2

for (i in 2:totalg) {
    # update precision
    kpost[i,] <- rgamma(2, postshape, hyperB + c(uRuHalf, vvHalf))

    # find eta tilde 
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

    # simulate proposal
    x_ <- c(rmvnorm.canonical(1, b, Q, Rstruct = struct))
    upost[i,] <- x_[1:n]; eta_ <- x_[1:n + n]; vpost[i,] <- eta_ - upost[i,]
    uRuHalf_ <- t(upost[i,]) %*% (R %*% upost[i,]) / 2
    vvHalf_ <- t(vpost[i,]) %*% vpost[i,] / 2

    # calculate acceptance probability
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

    # accept or reject proposal
    if (log(runif(1)) < logAlpha) {    
        uRuHalf <- uRuHalf_;  vvHalf <- vvHalf_
        eta <- eta_; b <- b_; C <- C_; accept[i] <- 1
    } else{                            
        accept[i] <- 0; upost[i,] <- upost[i-1,]; vpost[i,] <- vpost[i-1,]
    } 

    # progress report
    if(i%%500 == 0) cat(paste(i / 50, "%\n", sep = "" ))
}
 


# check values
head(eta)
tail(b)
head(C)
tail(accept)
tail(upost[30,])
tail(vpost[30,])
sum(accept[-1])
sum(upost)


if(FALSE){
## remove burnin 
kpost <- kpost[-seq(burnin), ]
upost <- upost[-seq(burnin), ]
vpost <- vpost[-seq(burnin), ]
accept <- accept[-seq(burnin)]

## thinning
index <- c(TRUE, rep(FALSE, (thin - 1)))
kpost <- kpost[index,]
upost <- upost[index,]
vpost <- vpost[index,]
accept <- accept[index]

## acceptance rate
mean(accept)
plot(accept, yaxt = "n")
axis(2, at = c(0,1), label = c("no", "yes"), las = 2)

## trace and mixing plots for the precision parameters
## kappa_u and kappa_v
grid.newpage()
grid_trace2(kpost, chain1_lab = expression(log~kappa[u]),
            chain2_lab = expression(log~kappa[v]), log = TRUE)

## summary statistics of 
## kappa_u and kappa_v
apply(kpost, 2, summary)           


par(mfrow = c(1,2))
## standardized mortality ratios
germany.plot(log(Y/E), main = "SMR")              
## estimated relative log-risk
germany.plot(apply(upost, 2, mean), main = "U | Y, hyper-priors") 

}

options( echo=TRUE)

