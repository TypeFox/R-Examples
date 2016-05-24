source("helper-diversitree.R")

context("Mkn")

## Simulate a tree and character distribution.  This is on a birth-death
## tree, with high rates of character evolution and an asymmetry in the
## character transition rates.
pars <- c(.1, .1, .03, .03, .1, .2)
set.seed(3)
phy <- trees(pars, "bisse", max.taxa=25, max.t=Inf, x0=0)[[1]]

## Maximum likelihood parameter estimation:
p <- c(.1, .1) # initial parameter guess
lik <- make.mk2(phy, phy$tip.state)
fit.mk2 <- find.mle(lik, p)
expect_that(fit.mk2$lnLik, equals(-10.90569519216388))

## Now, start the hard work...
control.p.d <- list(method="pij", backend="deSolve")
control.p.g <- list(method="pij", backend="gslode", compiled=FALSE)
control.p.G <- list(method="pij", backend="gslode", compiled=TRUE)
control.o.d <- list(method="ode", backend="deSolve")
control.o.g <- list(method="ode", backend="gslode", compiled=FALSE)
control.o.G <- list(method="ode", backend="gslode", compiled=TRUE)

states <- phy$tip.state + 1

lik.p.d <- make.mkn(phy, states, 2, control=control.p.d)
lik.p.g <- make.mkn(phy, states, 2, control=control.p.g)
lik.p.G <- make.mkn(phy, states, 2, control=control.p.G)
lik.o.d <- make.mkn(phy, states, 2, control=control.o.d)
lik.o.g <- make.mkn(phy, states, 2, control=control.o.g)
lik.o.G <- make.mkn(phy, states, 2, control=control.o.G)

ll <- fit.mk2$lnLik
expect_that(lik.p.d(coef(fit.mk2)), equals(ll))
expect_that(lik.p.g(coef(fit.mk2)), equals(ll))
expect_that(lik.p.G(coef(fit.mk2)), equals(ll))
expect_that(lik.o.d(coef(fit.mk2)), equals(ll))
expect_that(lik.o.g(coef(fit.mk2)), equals(ll))
expect_that(lik.o.G(coef(fit.mk2)), equals(ll))

## Check against ape.
model <- matrix(c(0, 2, 1, 0), 2)
fit.ape <- ace(phy$tip.state, phy, "discrete", model=model, ip=p)
expect_that(lik(fit.ape$rates, root=ROOT.GIVEN, root.p=c(1,1)),
            equals(fit.ape$loglik))

set.seed(3)
k <- 200
u <- 3.1
d <- 2.9
states <- sim.character(phy, c(k, u, d), model="meristic", x0=4)

control.m.d <- list(backend="deSolve")
control.m.g <- list(backend="gslode", compiled=FALSE)
control.m.G <- list(backend="gslode", compiled=TRUE)

lik.o <- make.mkn(phy, states, k, strict=FALSE,
                    control=list(method="ode"))
lik.m.0 <- make.mkn.meristic(phy, states, k)
lik.m.d <- make.mkn.meristic(phy, states, k, control=control.m.d)
lik.m.g <- make.mkn.meristic(phy, states, k, control=control.m.g)
lik.m.G <- make.mkn.meristic(phy, states, k, control=control.m.G)

pars <- diversitree:::mkn.meristic.Q(c(u, d), k)
pars.o <- pars[row(pars) != col(pars)]
pars.m <- c(d, u)

ll.o <- lik.o(pars.o)
expect_that(lik.m.0(pars.m), equals(ll.o))
expect_that(lik.m.d(pars.m), equals(ll.o, tolerance=1e-7))
expect_that(lik.m.g(pars.m), equals(ll.o))
expect_that(lik.m.G(pars.m), equals(ll.o))

