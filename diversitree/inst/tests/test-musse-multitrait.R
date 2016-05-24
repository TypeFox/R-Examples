source("helper-diversitree.R")

equals7 <- function(...)
  equals(..., tolerance=1e-7)

context("MuSSE (multitrait)")

tr <- musse.multitrait.translate(2)
pars <- c(.10, .15, .20, .25, # lambda 00, 10, 01, 11
          .03, .03, .03, .03, # mu 00, 10, 01, 11
          .05, .05, .0,       # q00.10, q00.01, q00.11
          .05, .0,  .05,      # q10.00, q10.01, q10.11
          .05, .0,  .05,      # q01.00, q01.10, q01.11
          .0,  .05, .05)      # q11.00, q11.10, q11.01
set.seed(2)
phy <- tree.musse(pars, 60, x0=1)

states <- expand.grid(A=0:1, B=0:1)[phy$tip.state,]
rownames(states) <- phy$tip.label

## Note that transition from the original MuSSE basis to this basis is
## only possible in general when depth=n.trait and allow.multistep=TRUE
## (as only this generates a square matrix that is invertible).
## However, when it is possible to express the set of parameters in the
## new basis (as it is above), this can be done through a pseudoinverse
## (here, a left inverse).
pars2 <- drop(solve(t(tr) %*% tr) %*% t(tr) %*% pars)

## Going from our new basis to the original MuSSE parameters is always
## straightforward.  This is done automatically in the likelihood
## function.
expect_that(c(tr %*% pars2), equals(pars, check.attributes=FALSE))

## This shows that the two traits act additively on speciation rate
## (lambdaAB is zero), that there is no effect of any trait on
## extinction (the only nonzero mu parameter is mu0) and transition
## rates for one trait are unaffected by other traits (the only nonzero
## q parameters are the qXij.0 parameters; qXij.Y parameters are all
## zero).

control.0 <- list()
control.d <- list(backend="deSolve")
control.g <- list(backend="gslode", compiled=FALSE)
control.G <- list(backend="gslode", compiled=TRUE)

## Here is our new MuSSE function parametrised as a multi-trait
## function:
lik.m.0 <- make.musse.multitrait(phy, states)
lik.m.d <- make.musse.multitrait(phy, states, control=control.d)
lik.m.g <- make.musse.multitrait(phy, states, control=control.g)
lik.m.G <- make.musse.multitrait(phy, states, control=control.G)

## Here "b" stands for base.
lik.b.0 <- make.musse(phy, phy$tip.state, 4)
lik.b.d <- make.musse(phy, phy$tip.state, 4, control=control.d)
lik.b.g <- make.musse(phy, phy$tip.state, 4, control=control.g)
lik.b.G <- make.musse(phy, phy$tip.state, 4, control=control.G)

liks <- list(lik.m.0, lik.m.d, lik.m.g, lik.m.G)

for ( f in liks ) {
  expect_that(f(pars2),
              equals7(lik.b.0(pars)))
  expect_that(f(pars2, root=ROOT.FLAT),
              equals7(lik.b.0(pars,  root=ROOT.FLAT)))
  expect_that(f(pars2, root=ROOT.OBS),
              equals7(lik.b.0(pars,  root=ROOT.OBS)))
  expect_that(f(pars2, root=ROOT.FLAT, condition.surv=TRUE),
              equals7(lik.b.0(pars,  root=ROOT.FLAT, condition.surv=TRUE)))
  expect_that(f(pars2, root=ROOT.OBS,  condition.surv=TRUE),
              equals7(lik.b.0(pars,  root=ROOT.OBS,  condition.surv=TRUE)))
}
