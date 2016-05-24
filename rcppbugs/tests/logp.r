
mu <- 0
tau <- 0.001
b <- mcmc.normal(rnorm(10),mu=mu,tau=tau)
all.equal(logp(b),
          sum(dnorm(as.vector(b),mean=mu,sd=1/sqrt(tau),log=TRUE)),
          tolerance=0.01,
          check.attributes = FALSE)

