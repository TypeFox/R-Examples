source("helper-diversitree.R")

context("BiSSE (split)")

## Different control parameters
control.d <- list(backend="deSolve")
control.g <- list(backend="gslode", compiled=FALSE)
control.G <- list(backend="gslode", compiled=TRUE)
control.i <- list(backend="invalid_backend")

## First, simulate the tree:
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
set.seed(546)
phy <- tree.bisse(pars, max.taxa=30, x0=0)

lik.b <- make.bisse(phy, phy$tip.state)
lik.b(pars) # -93.62479
expect_that(lik.b(pars), equals(-93.6247932708028))

nodes <- c("nd15", "nd18", "nd26")
nodes.i <- match(nodes, phy$node.label) + length(phy$tip.label)
split.t <- Inf # TODO: Drop
lik.s.0 <- make.bisse.split(phy, phy$tip.state, nodes.i, split.t)
lik.s.d <- make.bisse.split(phy, phy$tip.state, nodes.i, split.t,
                            control=control.d)
lik.s.g <- make.bisse.split(phy, phy$tip.state, nodes.i, split.t,
                            control=control.g)
lik.s.G <- make.bisse.split(phy, phy$tip.state, nodes.i, split.t,
                            control=control.G)

pars.s <- rep(pars, 4)

## Run the likelihod calculation:
expect_that(lik.s.0(pars.s), equals(lik.b(pars)))
expect_that(lik.s.d(pars.s), equals(lik.b(pars), tolerance=1e-7))
expect_that(lik.s.g(pars.s), equals(lik.b(pars)))
expect_that(lik.s.G(pars.s), equals(lik.b(pars)))

## Check the node identity...
lik.s2 <- make.bisse.split(phy, phy$tip.state, nodes, split.t)
expect_that(lik.s.0(pars.s), is_identical_to(lik.s2(pars.s)))

set.seed(1)
pars2 <- pars + runif(6, 0, .1)
pars2.s <- rep(pars2, 4)
expect_that(lik.b(pars2), equals(-94.162724105832))
expect_that(lik.s.0(pars2.s), equals(lik.b(pars2)))
expect_that(lik.s.d(pars2.s), equals(lik.b(pars2), tolerance=1e-7))
expect_that(lik.s.g(pars2.s), equals(lik.b(pars2)))
expect_that(lik.s.G(pars2.s), equals(lik.b(pars2)))

pars3.s <- pars + runif(length(pars.s), 0, .1)
ll <- -94.18230685862810
expect_that(lik.s.0(pars3.s), equals(ll))
expect_that(lik.s.d(pars3.s), equals(ll))
expect_that(lik.s.g(pars3.s), equals(ll))
expect_that(lik.s.G(pars3.s), equals(ll))

unresolved <-
  data.frame(tip.label=c("sp12", "sp32", "sp9", "sp22", "sp11"),
             Nc=c(2,5,3,2,5), n0=c(1, 4, 3, 2, 4), n1=c(1, 1, 0, 0, 1))

## Plain BiSSE with unresolved clades:
lik.u.b.0 <- make.bisse(phy, phy$tip.state, unresolved=unresolved)
lik.u.b.d <- make.bisse(phy, phy$tip.state, unresolved=unresolved,
                      control=control.d)
lik.u.b.g <- make.bisse(phy, phy$tip.state, unresolved=unresolved,
                      control=control.g)
lik.u.b.G <- make.bisse(phy, phy$tip.state, unresolved=unresolved,
                      control=control.G)

ll <- -139.36879978895746
expect_that(lik.u.b.0(pars), equals(ll))
expect_that(lik.u.b.d(pars), equals(ll))
expect_that(lik.u.b.g(pars), equals(ll))
expect_that(lik.u.b.G(pars), equals(ll))

## Split BiSSE with unresolved clades:
lik.u.s.0 <- make.bisse.split(phy, phy$tip.state, nodes, split.t,
                              unresolved=unresolved)
lik.u.s.d <- make.bisse.split(phy, phy$tip.state, nodes, split.t,
                              unresolved=unresolved,
                              control=control.d)
lik.u.s.g <- make.bisse.split(phy, phy$tip.state, nodes, split.t,
                              unresolved=unresolved,
                              control=control.g)
lik.u.s.G <- make.bisse.split(phy, phy$tip.state, nodes, split.t,
                              unresolved=unresolved,
                              control=control.G)

expect_that(lik.u.s.0(pars.s), equals(lik.u.b.0(pars)))
expect_that(lik.u.s.d(pars.s), equals(lik.u.b.0(pars)))
expect_that(lik.u.s.g(pars.s), equals(lik.u.b.0(pars)))
expect_that(lik.u.s.G(pars.s), equals(lik.u.b.0(pars)))

expect_that(lik.u.s.0(pars2.s), equals(lik.u.b.0(pars2)))
expect_that(lik.u.s.d(pars2.s), equals(lik.u.b.0(pars2), tolerance=1e-7))
expect_that(lik.u.s.g(pars2.s), equals(lik.u.b.0(pars2)))
expect_that(lik.u.s.G(pars2.s), equals(lik.u.b.0(pars2)))

ll <- -127.78056192052314
expect_that(lik.u.s.0(pars3.s), equals(ll))
expect_that(lik.u.s.d(pars3.s), equals(ll, tolerance=1e-7))
expect_that(lik.u.s.g(pars3.s), equals(ll))
expect_that(lik.u.s.G(pars3.s), equals(ll))
