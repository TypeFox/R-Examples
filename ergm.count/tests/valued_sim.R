#  File tests/valued_sim.R in package ergm.count, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2008-2016 Statnet Commons
#######################################################################
library(ergm.count)
load("testnet3d.RData")

set.seed(0)

## Poisson-reference
cat("======== Poisson-reference ERGM with mean 2, TNT proposal\n")
s <- simulate(testnet3d~sum, nsim=1000, reference=~Poisson, response="w", coef=log(2), statsonly=TRUE, control=control.simulate.formula(MCMC.prop.args=list(p0=0.8)))
test <- approx.hotelling.diff.test(s/6,mu0=2)
if(test$p.value<0.001) {print(test); stop("Simulation test failed.")}

cat("======== Poisson-reference ERGM with mean 2, random proposal\n")
s <- simulate(testnet3d~sum, nsim=1000, reference=~Poisson, response="w", coef=log(2), statsonly=TRUE, control=control.simulate.formula(MCMC.prop.weights="random"))
test <- approx.hotelling.diff.test(s/6,mu0=2)
if(test$p.value<0.001) {print(test); stop("Simulation test failed.")}

cat("======== Poisson-reference ERGM with mean 2, zero-modified proposal\n")
s <- simulate(testnet3d~sum, nsim=1000, reference=~Poisson, response="w", coef=log(2), statsonly=TRUE, control=control.simulate.formula(MCMC.prop.weights="0inflated", MCMC.prop.args=list(p0=0.8)))
test <- approx.hotelling.diff.test(s/6,mu0=2)
if(test$p.value<0.001) {print(test); stop("Simulation test failed.")}

## Geometric-reference
cat("======== Geometric-reference ERGM with mean 2\n")
s <- simulate(testnet3d~sum, nsim=1000, reference=~Geometric, response="w", coef=log(2/3), statsonly=TRUE)
test <- approx.hotelling.diff.test(s/6,mu0=2)
if(test$p.value<0.001) {print(test); stop("Simulation test failed.")}

## Binomial-reference
cat("======== Binomial-reference ERGM with mean 5 trials and probability of success 0.4 for a mean 2\n")
s <- simulate(testnet3d~sum, nsim=1000, reference=~Binomial(5), response="w", coef=log(.4/.6), statsonly=TRUE)
test <- approx.hotelling.diff.test(s/6,mu0=2)
if(test$p.value<0.001) {print(test); stop("Simulation test failed.")}
