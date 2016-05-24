source("helper-diversitree.R")

## Will be useful to have a function for testing tolerance to within
## 1e-7, as that works out to be how accurate most things actually
## are.
equals7 <- function(...)
  equals(..., tolerance=1e-7)

context("MuSSE")

## 1: BiSSE equivalence
pars <- c(.1, .2, .03, .04, 0.05, 0.1)
set.seed(2)
phy <- tree.musse(pars, 20, x0=1)

## Show that the likelihood functions give the same answers.
lik.b <- make.bisse(phy, phy$tip.state-1)
lik.m <- make.musse(phy, phy$tip.state, 2)

## Notice that default argument names are different between BiSSE and
## MuSSE, but that the order is the same.
argnames(lik.b) # BiSSE: 0/1
argnames(lik.m) # MuSSE: 1/2

## 2: A 3-state example where movement is only allowed between
## neighbouring states (1 <-> 2 <-> 3), and where speciation and
## extinction rates increase moving from 1 -> 2 -> 3:

## You can get the expected argument order for any number of states
## this way (sorry - clunky).  The help file also lists the order.
diversitree:::default.argnames.musse(3)

## Here are the parameters:
pars <- c(.1,  .15,  .2,  # lambda 1, 2, 3
          .03, .045, .06, # mu 1, 2, 3
          .05, 0,         # q12, q13
          .05, .05,       # q21, q23
          0,   .05)       # q31, q32

set.seed(2)
phy <- tree.musse(pars, 30, x0=1)

## The states are numbered 1:3, rather than 0:1 in bisse.
states <- phy$tip.state
table(states)

## 2: Likelihood
## Making a likelihood function is basically identical to bisse.  The
## third argument needs to be the number of states.  In a future
## version this will probably be max(states), but there are some
## pitfalls about this that I am still worried about.
lik <- make.musse(phy, states, 3)

## Here are the arguments.  Even with three states, this is getting
## ridiculous.
argnames(lik)

lik.d <- make.musse(phy, states, 3, control=list(backend="deSolve"))
lik.g <- make.musse(phy, states, 3, control=list(compiled=FALSE))
set.seed(1)
p <- runif(length(pars), max=.05)

expect_that(lik.d(p), equals(lik(p)))
expect_that(lik.g(p), equals(lik(p)))

## Start with a fully constrained model, but still enforcing stepwise
## changes (disallowing 1 <-> 3 shifts)
lik.base <- constrain(lik, lambda2 ~ lambda1, lambda3 ~ lambda1,
                      mu2 ~ mu1, mu3 ~ mu1,
                      q13 ~ 0, q21 ~ q12, q23 ~ q12, q31 ~ 0, q32 ~ q12)

p <- starting.point.musse(phy, 3)

fit.base <- find.mle(lik.base, p[argnames(lik.base)])
expect_that(fit.base$lnLik, equals(-110.634755712668))

## Now, allow the speciation rates to vary:
lik.lambda <- constrain(lik, mu2 ~ mu1, mu3 ~ mu1,
                        q13 ~ 0, q21 ~ q12, q23 ~ q12, q31 ~ 0, q32 ~ q12)
fit.lambda <- find.mle(lik.lambda, p[argnames(lik.lambda)])
expect_that(fit.lambda$lnLik, equals(-110.211741798546))
