source("helper-diversitree.R")

## Will be useful to have a function for testing tolerance to within
## 1e-7, as that works out to be how accurate most things actually
## are.
equals7 <- function(...)
  equals(..., tolerance=1e-7)

context("MuSSE (time-dependent)")

pars <- c(.1,  .15,  .2,  # lambda 1, 2, 3
          .03, .045, .06, # mu 1, 2, 3
          .05, 0,         # q12, q13
          .05, .05,       # q21, q23
          0,   .05)       # q31, q32
set.seed(2)
phy <- tree.musse(pars, 30, x0=1)

control.0 <- list()
control.d <- list(backend="deSolve")
control.g <- list(backend="gslode", compiled=FALSE)
control.G <- list(backend="gslode", compiled=TRUE)

lik.m <- make.musse(phy, phy$tip.state, 3)
functions <- rep(c("sigmoid.t", "constant.t"), c(3, 9))

lik.0 <- make.musse.t(phy, phy$tip.state, 3, functions, control=control.0)
lik.d <- make.musse.t(phy, phy$tip.state, 3, functions, control=control.d)
lik.g <- make.musse.t(phy, phy$tip.state, 3, functions, control=control.g)
lik.G <- make.musse.t(phy, phy$tip.state, 3, functions, control=control.G)

## First, evaluate the functions with no time effect and check that they
## are the same as the base MuSSE model
t <- max(branching.times(phy))/2
p.t <- c(pars[1], pars[1], t, .5,
         pars[2], pars[2], t, .5,
         pars[3], pars[3], t, .5,
         pars[4:12])

ll <- lik.m(pars)
expect_that(ll, equals(-110.836441114592))

expect_that(lik.0(p.t), equals(ll))
expect_that(lik.d(p.t), equals7(ll))
expect_that(lik.g(p.t), equals(ll))
expect_that(lik.G(p.t), equals(ll))

## With some random parameters, still no time dependence:
set.seed(1)
pars2 <- pars + runif(length(pars), 0, .2)
p2.t <- c(pars2[1], pars2[1], t, .5,
          pars2[2], pars2[2], t, .5,
          pars2[3], pars2[3], t, .5,
          pars2[4:12])

ll <- lik.m(pars2)
expect_that(ll, equals(-118.952881268597))
## And expect that the time dependent versions agree again:
expect_that(lik.0(p2.t), equals(ll))
expect_that(lik.d(p2.t), equals7(ll))
expect_that(lik.g(p2.t), equals(ll))
expect_that(lik.G(p2.t), equals(ll))

## Add some time dependence:
p3.t <- p2.t
p3.t[2] <- p3.t[2] + .2
ll <- lik.0(p3.t)
expect_that(ll, equals(-118.752045964306))
expect_that(lik.d(p3.t), equals7(ll))
expect_that(lik.g(p3.t), equals(ll))
expect_that(lik.G(p3.t), equals(ll))

## TODO: Repeat but for variation in Q matrix?
