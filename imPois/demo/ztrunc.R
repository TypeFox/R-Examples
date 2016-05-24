#
# This demo shows comparison between numerical methods on estimating
# parameter of interest when zero-truncated Poisson sampling model 
# is considered. 
# 
# updated on 2015.11.27

lambda <- 1
n <- 1e2
y <- rpois(n=n, lambda=lambda)
y <- y[y>0]

opaq <- evfn.ztrunc(y=y, pars=c(xi2=0, xi1=1, xi0=1))
opmh <- mh.ztrunc(y=y, pars=c(xi2=0, xi1=1, xi0=1))
# hist(op$chain)
opla <- lapprox.ztrunc(y=y, pars=c(xi2=0, xi1=1, xi0=1), const=1e3)

theta <- c(AQ=opaq$value, MH=opmh$value, LA=opla$value)
mu <- exp(theta)
cbind(theta, mu)

