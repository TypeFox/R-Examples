# portOpt3.R -- version 2010-12-11
# tangency portfolio
require(quadprog)

# create random returns
na  <- 20L  # number of assets 
ns  <- 60L  # number of observations
R   <- array(rnorm(ns*na, mean = 0.005, sd = 0.015), 
                  dim = c(ns,na))
mu  <- colMeans(R) # means
rf  <- 0.0001      # riskfree rate (about 2.5% pa)
mu2 <- mu - rf     # excess means

# set up matrices
Q <- cov(R)               # covariance matrix
B <- array(mu2, dim = c(1L,na))
b <- 1
result <- solve.QP(Dmat = Q,
                   dvec = rep(0,na),
                   Amat = t(B),
                   bvec = b,
                   meq  = 1L)

# rescale variables to obtain weights
w <- as.matrix(result$solution/sum(result$solution))

# compute sharpe ratio
SR <- t(w) %*% mu2 / sqrt(t(w) %*% Q %*% w)
sum(w)             # check budget constraint
t(w) %*% mu >= rf  # check return constraint

# test 1: regression approach from Britten-Jones (1999)
R2   <- R - rf
ones <- array(1, dim = c(ns,1L))
solR <- lm(ones~-1 + R2)
w2 <- coef(solR); w2 <- w2/sum(w2)
# ... w2 should be the same as w
stopifnot(all.equal(as.numeric(w),as.numeric(w2)))

# test 2: no inequality constraints >> solve FOC 
w3 <- solve(Q,mu2)
w3 <- w3/sum(w3)
# ... w3 should be the same as w2 and w
stopifnot(all.equal(as.numeric(w2),as.numeric(w3)))