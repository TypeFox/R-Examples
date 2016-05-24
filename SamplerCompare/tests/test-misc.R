# From SamplerCompare, (c) 2010 Madeleine Thompson

# This script contains a couple miscellaneous tests.

library(SamplerCompare)

# Make sure the R-C glue for distributions works by creating C and
# R versions of the same 2-D Gaussian and making sure their log density
# and its gradient agree on a randomly chosen point.

N2.C <- make.c.dist(2, 'Gauss2-C', 'Gauss2_log_dens', c(1, 2, 0.8), mean=c(1,2))
N2.R <- make.gaussian(c(1,2), rho=0.8)
x <- runif(2)

stopifnot(abs(N2.C$log.density(x)-N2.R$log.density(x)) < 1e-08)
stopifnot(max(abs(N2.C$grad.log.density(x)-N2.R$grad.log.density(x))) < 1e-08)

# Make sure multivariate Gamma distributions have a log density
# that matches both dgamma and their gradient.

k <- c(1,2)
theta <- c(3,4)
ds <- make.mv.gamma.dist(k, theta)
x <- runif(2)
stopifnot(abs(sum(dgamma(x, k, scale=theta, log=TRUE)) -
              ds$log.density(x)) < 1e-5)
check.dist.gradient(ds, x)
