source("helper-diversitree.R")

## Will be useful to have a function for testing tolerance to within
## 1e-7, as that works out to be how accurate most things actually
## are.
equals7 <- function(...)
  equals(..., tolerance=1e-7)

context("Birth-death")

set.seed(1)
phy <- trees(c(.1, .03), "bd", max.taxa=25)[[1]]

control.n <- list(method="nee")
control.d <- list(method="ode", backend="deSolve")
control.g <- list(method="ode", backend="gslode", compiled=FALSE)
control.G <- list(method="ode", backend="gslode", compiled=TRUE)

lik.0 <- make.bd(phy)
lik.n <- make.bd(phy, control=control.n)
lik.d <- make.bd(phy, control=control.d)
lik.g <- make.bd(phy, control=control.g)
lik.G <- make.bd(phy, control=control.G)

p <- c(.1, .03)
ll <- -22.50266563415211
expect_that(lik.0(p), equals(ll))
expect_that(lik.n(p), equals(ll))
expect_that(lik.d(p), equals(ll, tolerance=1e-6))
expect_that(lik.g(p), equals(ll))
expect_that(lik.G(p), equals(ll))

set.seed(1)
p <- c(.1, .03) + runif(2, 0, .2)
ll <- -23.59180370809460
expect_that(lik.0(p), equals(ll))
expect_that(lik.n(p), equals(ll))
expect_that(lik.d(p), equals(ll, tolerance=1e-6))
expect_that(lik.g(p), equals(ll))
expect_that(lik.G(p), equals(ll))

ll <- -27.26771051981412
expect_that(lik.0(p, condition.surv=FALSE), equals(ll))
expect_that(lik.n(p, condition.surv=FALSE), equals(ll))
expect_that(lik.d(p, condition.surv=FALSE), equals(ll, tolerance=1e-7))
expect_that(lik.g(p, condition.surv=FALSE), equals(ll))
expect_that(lik.G(p, condition.surv=FALSE), equals(ll))

fit.ape <- birthdeath(phy)
p.ape <- rev(fit.ape$para) # zero extinction..

fit.dt <- find.mle(lik.0, p, method="optim", lower=0)
expect_that(coef(fit.dt), equals(p.ape, tolerance=1e-5,
                                 check.attributes=FALSE))
expect_that(fit.dt$lnLik, equals(-fit.ape$dev/2))

lik.s.0 <- make.bd(phy, sampling.f=.5)
lik.s.n <- make.bd(phy, sampling.f=.5, control=control.n)
lik.s.d <- make.bd(phy, sampling.f=.5, control=control.d)
lik.s.g <- make.bd(phy, sampling.f=.5, control=control.g)
lik.s.G <- make.bd(phy, sampling.f=.5, control=control.G)

ll <- -23.81317463410290
expect_that(lik.s.0(p), equals(ll))
expect_that(lik.s.n(p), equals(ll))
expect_that(lik.s.d(p), equals(ll, tolerance=1e-6))
expect_that(lik.s.g(p), equals(ll))
expect_that(lik.s.G(p), equals(ll))

ll <- -27.7338143635087
expect_that(lik.s.0(p, condition.surv=FALSE), equals(ll))
expect_that(lik.s.n(p, condition.surv=FALSE), equals(ll))
expect_that(lik.s.d(p, condition.surv=FALSE), equals(ll, tolerance=1e-6))
expect_that(lik.s.g(p, condition.surv=FALSE), equals(ll))
expect_that(lik.s.G(p, condition.surv=FALSE), equals(ll))

