source("helper-diversitree.R")

## Will be useful to have a function for testing tolerance to within
## 1e-7, as that works out to be how accurate most things actually
## are.
equals7 <- function(...)
  equals(..., tolerance=1e-7)

context("BiSSE (time-dependent)")

set.seed(4)
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
phy <- tree.bisse(pars, max.t=30, x0=0)

control.0 <- list()
control.d <- list(backend="deSolve")
control.g <- list(backend="gslode", compiled=FALSE)
control.G <- list(backend="gslode", compiled=TRUE)

lik.b <- make.bisse(phy, phy$tip.state)
functions <- rep(c("sigmoid.t", "constant.t"), c(2, 4))

lik.0 <- make.bisse.t(phy, phy$tip.state, functions, control=control.0)
lik.d <- make.bisse.t(phy, phy$tip.state, functions, control=control.d)
lik.g <- make.bisse.t(phy, phy$tip.state, functions, control=control.g)
lik.G <- make.bisse.t(phy, phy$tip.state, functions, control=control.G)

## First, evaluate the functions with no time effect and check that they
## are the same as the base BiSSE model
t <- max(branching.times(phy))/2
p.t <- c(pars[1], pars[1], t, .5,
         pars[2], pars[2], t, .5,
         pars[3:6])

ll <- lik.b(pars)
expect_that(ll, equals(-159.709959436709))

expect_that(lik.0(p.t), equals(ll))
expect_that(lik.d(p.t), equals7(ll))
expect_that(lik.g(p.t), equals(ll))
expect_that(lik.G(p.t), equals(ll))

## With some random parameters, still no time dependence:
set.seed(1)
pars2 <- pars + runif(6, 0, .2)
p2.t <- c(pars2[1], pars2[1], t, .5,
          pars2[2], pars2[2], t, .5,
          pars2[3:6])

ll <- lik.b(pars2)
expect_that(ll, equals(-189.608098958571))
## And expect that the time dependent versions agree again:
expect_that(lik.0(p2.t), equals(ll))
expect_that(lik.d(p2.t), equals7(ll))
expect_that(lik.g(p2.t), equals(ll))
expect_that(lik.G(p2.t), equals(ll))

## Add some time dependence:
p3.t <- p2.t
p3.t[2] <- p3.t[2] + .2
ll <- lik.0(p3.t)
expect_that(ll, equals(-190.581103942523))
expect_that(lik.d(p3.t), equals7(ll))
expect_that(lik.g(p3.t), equals(ll))
expect_that(lik.G(p3.t), equals(ll))

