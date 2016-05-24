source("helper-diversitree.R")

## Will be useful to have a function for testing tolerance to within
## 1e-7, as that works out to be how accurate most things actually
## are.
equals7 <- function(...)
  equals(..., tolerance=1e-7)

context("BiSSE")

## Basic simulated tree, 30 taxa.
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
set.seed(4)
phy <- tree.bisse(pars, max.t=30, x0=0)

## Different control parameters
control.d <- list(backend="deSolve")
control.g <- list(backend="gslode", compiled=FALSE)
control.G <- list(backend="gslode", compiled=TRUE)
control.i <- list(backend="invalid_backend")

lik.0 <- make.bisse(phy, phy$tip.state)
lik.d <- make.bisse(phy, phy$tip.state, control=control.d)
lik.g <- make.bisse(phy, phy$tip.state, control=control.g)
lik.G <- make.bisse(phy, phy$tip.state, control=control.G)

## Check that invalid backends trigger error:
expect_that(make.bisse(phy, phy$tip.state, control=control.i),
            throws_error())

ll1 <- -159.709955958896
expect_that(lik.0(pars), equals7(ll1))
expect_that(lik.d(pars), equals7(ll1))
expect_that(lik.g(pars), equals7(ll1))
expect_that(lik.G(pars), equals7(ll1))

set.seed(1)
pars2 <- pars + runif(6, 0, .1)

opts <- data.frame(surv=rep(c(TRUE, FALSE), each=3),
                   root=rep(c(ROOT.FLAT, ROOT.OBS, ROOT.EQUI), 2),
                   ll=c(-173.403159254299,
                     -173.083837291932,
                     -173.57484784988,
                     -176.564186725815, 
                     -175.939144822762,
                     -176.840650405926))

liks <- list(lik.0, lik.d, lik.g, lik.G)

## To generate the list above:
## dput(mapply(function(r, s) lik.0(pars2, root=r, condition.surv=s),
##             opts$root, opts$surv))
for ( f in liks )
  for ( i in seq_len(nrow(opts)) )
    expect_that(f(pars2, root=opts$root[i], condition.surv=opts$surv[i]),
                equals7(opts$ll[i]))

## Here, generate an unresolved description that is the same as the
## individual tips (i.e., single species, in observed state)
unresolved <- data.frame(tip.label=phy$tip.label, Nc=1,
                         n0=1-phy$tip.state, n1=phy$tip.state)
states <- phy$tip.state
states[] <- NA

lik.u.0 <- make.bisse(phy, states, unresolved=unresolved)
lik.u.d <- make.bisse(phy, states, unresolved=unresolved,
                      control=control.d)
lik.u.g <- make.bisse(phy, states, unresolved=unresolved,
                      control=control.g)
lik.u.G <- make.bisse(phy, states, unresolved=unresolved,
                      control=control.G)

expect_that(lik.u.0(pars), equals7(lik.0(pars)))
expect_that(lik.u.d(pars), equals7(lik.0(pars)))
expect_that(lik.u.g(pars), equals7(lik.0(pars)))
expect_that(lik.u.G(pars), equals7(lik.0(pars)))

expect_that(lik.u.0(pars2), equals(lik.0(pars2), tolerance=3e-6))
expect_that(lik.u.d(pars2), equals(lik.0(pars2), tolerance=3e-6))
expect_that(lik.u.g(pars2), equals(lik.0(pars2), tolerance=3e-6))
expect_that(lik.u.G(pars2), equals(lik.0(pars2), tolerance=3e-6))
