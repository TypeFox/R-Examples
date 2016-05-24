# Compare naive and uniform random sampling from beta + eta < 1
# Uniform sampling is implemented in blim(..., randinit = TRUE).
#
# Last mod: Sep/12/2013, FW

nsamples <- 12000

## The following naive algorithm samples non-uniformly from beta + eta < 1

beta <- runif(nsamples)
 eta <- runif(nsamples, max=1 - beta)

plot(eta ~ beta, pch=".", main="Naive sampling")

## The following algorithm samples uniformly from beta + eta < 1

beta <- runif(nsamples)
 eta <- runif(nsamples)
beta <- ifelse(beta + eta < 1, beta, 1 - beta)
 eta <- ifelse(beta + eta < 1,  eta, 1 -  eta)

plot(eta ~ beta, pch=".", main="Uniform sampling")

