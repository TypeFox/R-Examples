# This is file ../spam/tests/demo_jss15-Leroux.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     

options( echo=FALSE)
library( spam, warn.conflict=FALSE)




# This is the MCMC sampler presented in Section 3.1 of the article:
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
# This takes only 1-2 minutes of computation time.
#
# In the JSS article we used:
# number of generated samples: 300'000
# number of burnin samples:     15'000
# thinning:                         20

## invisible(readline(prompt = "Type  <Return>\t to continue : "))

# SETUP:
spam.options(structurebased=TRUE)
library("truncdist")


# Lereux model - one block, no intercept

# load data
data(Oral)
E <- Oral$E
Y <- Oral$Y
n <- 544
A <- adjacency.landkreis(system.file("demodata/germany.adjacency", package="spam"))
A <- as.matrix(A)

set.seed(2)

# tuning parameters
# (influence acceptance rate)
lambda.proposal.sd <- 0.070992


# sampler length, burnin and thinning
totaln <- 100
burnin <- 500
thin <- 10

# variable to store samples
upost <- array(NA, c(totaln, n))
kpost <- rep(NA,totaln)
lpost <- rep(NA,totaln)
accept <- array(0, c(totaln, 2),
                list(NULL, c("U", "lambda")))

# initial values
upost[1,] <- rep(c(.1,-.1), 544/2)
kpost[1] <- 15
lpost[1] <- .9
accept[1,] <- 1

# precaluclate some quantities
R <- precmat.IGMRFirreglat(A)
eigenR <- eigen(R); eigenR.value <- eigenR$values; 

Q <- (1 - lpost[1]) * diag.spam(544) + lpost[1] * R
Q.det <- sum(log(lpost[1]* eigenR.value + 1 - lpost[1]))
Q.struct <- chol.spam(Q)

postshape <- 0.5 * n - 1

# start sampler
for (i in 2:totaln) {

    ## u
    ## find tilde u
    u.tilde <- upost[i-1,]
    C <- E * exp(u.tilde)
    B <- Y + (u.tilde - 1) * C
    Q.tmp <- diag.spam(C) + kpost[i-1] * Q 
    u.tilde <- c(solve.spam(Q.tmp, B))
    
    C.tilde <- E * exp(u.tilde)
    B.tilde <- Y + (u.tilde - 1) * C.tilde
    Q.tilde <- diag.spam(C.tilde) + kpost[i-1] * Q
    
    u_ <- c(rmvnorm.canonical(1, B.tilde, Q.tilde, Rstruct = Q.struct))

    u.tilde_ <- u_
    C_ <- E * exp(u.tilde_)
    B_ <- Y + (u.tilde_- 1) * C_
    Q.tmp_ <- diag.spam(C_) + kpost[i-1] * Q 
    u.tilde_ <- c(solve.spam(Q.tmp_, B_))
  
    C.tilde_ <- E * exp(u.tilde_)
    B.tilde_ <- Y + (u.tilde_- 1) * C.tilde_
 
    log.alpha.u <- sum(Y * u_) - sum(E*exp(u_))  -
      sum(Y * upost[i-1,]) + sum(E*exp(upost[i-1,]) ) +
        sum(upost[i-1,] * B.tilde_) -
          .5 * t(upost[i-1,]) %*% (diag(C.tilde_) %*% upost[i-1,]) -
            sum(u_* B.tilde) + .5 * t(u_) %*% (diag(C.tilde) %*% u_) 
    
    if(exp(log.alpha.u) > runif(1)){
      upost[i,] <- u_
      accept[i,1] <- 1
    } else{
      upost[i,] <- upost[i-1,]
    }
    
    ## kappa
    kpost[i] <- rgamma(1, shape = postshape,
                       rate = .5 * upost[i,] %*% (Q %*% upost[i,]))
    
    ## lambda
    lambda_ <- rtrunc(n=1, spec="norm", a=0, b=1,  mean=lpost[i-1], sd=lambda.proposal.sd)
    Q_ <- (1 - lambda_) * diag.spam(544) + lambda_ * R
    Q.det_ <- sum(log(lambda_* eigenR.value + 1 - lambda_))

  
    alpha.lambda <- exp(.5 * (Q.det_ -
                              kpost[i]* upost[i,] %*% (Q_ %*% upost[i,]) -
                              Q.det +
                              kpost[i]* upost[i,] %*% (Q %*% upost[i,]) ))
                         
    if(alpha.lambda > runif(1)){
      lpost[i] <- lambda_
      Q <- Q_
      Q.det <- Q.det_
      accept[i,2] <- 1
    }else{
      lpost[i] <- lpost[i-1]
    }

    # progress report
    if(i%%500 == 0) cat(paste(i / 50, "%\n", sep = "" ))   

}

tail(lpost)
tail(kpost)
upost[i,100:150]


if(FALSE){
# remove burnin 
kpost <- kpost[-seq(burnin)]
upost <- upost[-seq(burnin), ]
lpost <- lpost[-seq(burnin)]
accept <- accept[-seq(burnin),]

# thinning
index <- c(TRUE, rep(FALSE, (thin - 1)))
kpost <- kpost[index]
upost <- upost[index,]
lpost <- lpost[index]
accept <- accept[index,]

# acceptance rate
apply(accept, 2, mean)
par(mfrow = c(1,2))
plot(accept[,1], yaxt = "n", main = expression(beta))
axis(2, at = c(0,1), label = c("no", "yes"), las = 2)
plot(accept[,2], yaxt = "n", main = expression(lambda))
axis(2, at = c(0,1), label = c("no", "yes"), las = 2)

# trace and mixing plots for the precision parameters
# kappa and lambda
grid.newpage();
grid_trace2(log(kpost), lpost,
            chain1_lab = expression(log~kappa),
            chain2_lab = expression(lambda))

# summary statistics of 
# kappa and lambda
summary(kpost)           
summary(lpost)           


par(mfrow = c(1,2))
## standardized mortality ratios
germany.plot(log(Y/E), main = "SMR")              
## estimated relative log-risk
germany.plot(apply(upost, 2, mean), main = "U | Y, hyper-priors") 
}



options( echo=TRUE)
