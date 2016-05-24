# portReg.R -- version 2011-01-06
# create artifical data ('daily returns')
n  <- 100L  # number of observations
p  <- 10L   # number of assets
X  <- array(rnorm(n * p, mean = 0.001, sd = 0.01), 
           dim = c(n, p))
rf <- 0.0001  # riskfree rate (2.5% pa)
m  <- apply(X, 2, mean)  # means
m2   <- m - rf           # excess means

## (1) solve the problem with qp
require(quadprog)
aMat  <- as.matrix(m2); bVec  <- 1
zeros <- array(0, dim = c(p, 1L))
solQP <- solve.QP(cov(X), zeros, aMat, bVec, meq = 1L)
# rescale variables to obtain weights
w     <- solQP$solution/sum(solQP$solution)
# compute sharpe ratio
SR    <- t(w) %*% m2 / sqrt(t(w) %*% cov(X) %*% w)

## (2) solve with regression
X2     <- X - rf # excess returns
ones  <- array(1, dim = c(n, 1L))
# run regression
solR  <- lm(ones~-1 + X2)
# rescale variables to obtain weights
w2    <- coef(solR)
w2    <- w2/sum(w2)

## (3) solve first-order conditions
w3 <- solve(cov(X), m2)
# rescale
w3 <- w3/sum(w3)

## check they are the same
stopifnot(all.equal(as.vector(w),  as.vector(w2)))
stopifnot(all.equal(as.vector(w),  as.vector(w3)))
stopifnot(all.equal(as.vector(w2), as.vector(w3)))