# portOpt2.R -- version 2010-12-28
# mean--variance-efficient portfolio
require(quadprog)

# create random returns
na  <- 20L  # number of assets 
ns  <- 60L  # number of observations
R   <- array(rnorm(ns * na, mean = 0.005, sd = 0.015), 
             dim = c(ns, na))
m    <- colMeans(R)  # means
rd   <- mean(m)      # desired mean
wsup <- 0.1          # maximum holding size
winf <- 0.0          # minimum holding size

# set up matrices
Q <- 2 * cov(R)  
A <- array( 1, dim = c(1L, na))
a <- 1
B <- array(m, dim = c(1L, na))
B <- rbind(B,-diag(na),diag(na))
b <- rbind(rd, array(-wsup, dim = c(na,1L)),
               array( winf, dim = c(na,1L)))
result <- solve.QP(Dmat = Q,
                   dvec = rep(0,na),
                   Amat = t(rbind(A,B)),
                   bvec = rbind(a,b),
                   meq  = 1L)

w <- result$solution 
sum(w)         # check budget constraint
w %*% m >= rd  # check return constraint
summary(w)     # check holding size constraint