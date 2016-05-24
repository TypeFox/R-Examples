# portRegMV.R -- version 2011-01-06
# create artificial data with mean 0 and sd 5%
n <- 100L  # number of observations
p <- 10L   # number of assets
X <- array(rnorm(n * p, mean = 0, sd = 0.05), 
           dim = c(n, p))

# (1) solve with qp
require(quadprog)
aMat  <- array(1, dim = c(1L, p)); bVec  <- 1
zeros <- array(0, dim = c(p, 1L))
solQP <- solve.QP(cov(X), zeros, t(aMat), bVec, meq = 1L)
# ... and check solution
stopifnot(all.equal(as.numeric(var(X %*% solQP$solution)), 
          as.numeric(2 * solQP$value)))

# (2) regression
y  <- X[ ,1L]           # choose 1st asset as regressand
X2 <- X[ ,1L] - X[ ,2L:p] # choose 1st asset as regressand
solR <- lm(y ~ X2)
# compare results of regression with qp
# _weights from qp
as.vector(solQP$solution)
# _weights from regression
as.vector(c(1 - sum(coef(solR)[-1L]), coef(solR)[-1L]))
# variance of portfolio
stopifnot(all.equal(as.numeric(var(X %*% solQP$solution)), var(solR$residuals)))

# (3) solve first-order conditons
x <- solve(cov(X), numeric(p) + 1) # or any other constant != 0
# rescale
x <- x / sum(x)
stopifnot(all.equal(solQP$solution, x))