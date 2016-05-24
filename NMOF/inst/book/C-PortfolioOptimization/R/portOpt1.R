# portOpt1.R -- version 2010-12-11
# minimum-variance portfolio with budget constraint
require(quadprog)

# create random returns
na <- 10L  # number of assets 
ns <- 60L  # number of observations
R  <- array(rnorm(ns*na, mean = 0.005, sd = 0.015), 
            dim = c(ns,na))

# minimize variance
Q <- 2 * cov(R)
A <- rbind(rep(1,10L))
a <- 1
result <- solve.QP(Dmat = Q,
                   dvec = rep(0,10L),
                   Amat = t(A),
                   bvec = a,
                   meq  = 1L)

## check budget constraint and solution
w <- result$solution
sum(w)
stopifnot(all.equal(as.numeric(var(R %*% w)),result$value))