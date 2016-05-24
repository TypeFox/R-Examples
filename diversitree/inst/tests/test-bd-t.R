## TODO: Could do with more extensive tests here.
source("helper-diversitree.R")

## Will be useful to have a function for testing tolerance to within
## 1e-7, as that works out to be how accurate most things actually
## are.
equals7 <- function(...)
  equals(..., tolerance=1e-7)

context("Birth-death (time-dependent)")

set.seed(1)
pars <- c(.1, .03)
phy <- trees(pars, "bd", max.taxa=25)[[1]]

## Next, make three different likelihood functions: a "normal" one that
## uses the direct birth-death calculation, an "ode" based one (that
## uses numerical integration to compute the likelihood, and is
## therefore not exact), and one that is time-varying, but that the
## time-dependent functions are constant.t().
lik.nee <- make.bd(phy)
lik.ode <- make.bd(phy, control=list(method="ode"))
lik.t <- make.bd.t(phy, c("constant.t", "constant.t"))

expect_that(lik.nee(pars), equals(-22.50266563415211))

## ODE-based likelihood calculations are correct to about 1e-6.
expect_that(lik.ode(pars), equals(lik.nee(pars)))

## The ODE calculation agrees exactly with the time-varying (but
## constant) calculation.
expect_that(lik.ode(pars), is_identical_to(lik.t(pars)))

## Next, make a real case, where speciation is a linear function of
## time.
lik.t2 <- make.bd.t(phy, c("linear.t", "constant.t"))

## Confirm that this agrees with the previous calculations when the
## slope is zero
pars2 <- c(pars[1], 0, pars[2])
expect_that(lik.t2(pars2),  is_identical_to(lik.t(pars)))

set.seed(1)
pars3 <- pars2 + runif(length(pars2), 0, .2)
expect_that(lik.t2(pars3), equals(-105.437123043537))

## Time chunks
set.seed(1)
## BUG: tree.bd chokes with wrong length parameters...
phy <- tree.bd(c(.1, 0), max.taxa=30)

p0 <- c(.1, 0)
p1 <- c(.1, .05)
p0.t <- c(p0[1], 0, p0[2])
p1.t <- c(p1[1], 0, p1[2])
p2.t <- c(p1[1], .005, p1[2])

lik1 <- make.bd(phy, control=list(method="ode"))
ll0 <- lik1(p0)
ll1 <- lik1(p1)

lik2 <- make.bd.t(phy, rep("constant.t", 2))
expect_that(lik2(p0), equals(ll0))
expect_that(lik2(p1), equals(ll1))

## old version
ll2 <- -20.3329826871675
lik3 <- make.bd.t(phy, c("linear.t", "constant.t"))
expect_that(lik3(p0.t), equals(ll0))
expect_that(lik3(p1.t), equals(ll1))
expect_that(lik3(p2.t), equals(ll2))

## Also with deSolve:
lik4 <- make.bd.t(phy, c("linear.t", "constant.t"),
                  control=list(backend="deSolve"))
expect_that(lik4(p0.t), equals(ll0, tolerance=2e-7))
expect_that(lik4(p1.t), equals(ll1, tolerance=2e-7))
expect_that(lik4(p2.t), equals(ll2, tolerance=2e-7))

## Now, a spline fit:
t.max <- max(branching.times(phy)) * 1.001

## First, the full function (a sin wave down to the root of the
## tree).
sin.t <- function(t, y0, y1)
  y0 + (y1 - y0) * (sin(t/t.max*2*pi)+1)/2

p3.t <- c(p0[1], p0[1], p0[2])
p4.t <- c(p0[1], 2*p0[1], p0[2])

## Then, with a simple fit through this.
x <- seq(0, t.max, length.out=101)
y <- sin(x/t.max*2*pi)
spline.data <- list(t=x, y=y)

lik5 <- make.bd.t(phy, c("spline.t", "constant.t"),
                  spline.data=spline.data)
expect_that(lik5(p3.t), equals(ll0))
expect_that(lik5(p4.t), equals(-23.0593116707151))

