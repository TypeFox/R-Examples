source("helper-diversitree.R")

## Different control parameters
control.d <- list(backend="deSolve")
control.g <- list(backend="gslode", compiled=FALSE)
control.G <- list(backend="gslode", compiled=TRUE)
control.i <- list(backend="invalid_backend")

context("MuSSE (split)")

pars <- c(.1,  .15,  .2,  # lambda 1, 2, 3
          .03, .045, .06, # mu 1, 2, 3
          .05, 0,         # q12, q13
          .05, .05,       # q21, q23
          0,   .05)       # q31, q32

set.seed(2)
phy <- tree.musse(pars, 30, x0=1)

## Basic function for comparison:
lik.m <- make.musse(phy, phy$tip.state, 3)
ll <- -110.836441114592
expect_that(lik.m(pars), equals(ll))

## Split this phylogeny at three points: nd16 and nd25, splitting it
## into three chunks
nodes <- c("nd16", "nd25")

## To make a split XXSSE function, pass the node locations and times
## in.  Here, we'll use 'Inf' as the split time to mimick MEDUSA's
## behaviour of placing the split at the base of the branch subtended by
## a node.
lik.s.0 <- make.musse.split(phy, phy$tip.state, 3, nodes, split.t=Inf)
lik.s.d <- make.musse.split(phy, phy$tip.state, 3, nodes, split.t=Inf,
                            control=control.d)
lik.s.g <- make.musse.split(phy, phy$tip.state, 3, nodes, split.t=Inf,
                            control=control.g)
lik.s.G <- make.musse.split(phy, phy$tip.state, 3, nodes, split.t=Inf,
                            control=control.G)
                            
## The parameters must be a list of the same length as the number of
## partitions.  Partition '1' is the root partition, and partition i is
## the partition rooted at the node[i-1]:
## argnames(lik.s)

## Because we have two nodes, there are three sets of parameters.
## Replicate the original list to get a starting point for the analysis:
pars.s <- rep(pars, 3)
names(pars.s) <- argnames(lik.s.0)

ll <- lik.m(pars)
expect_that(lik.s.0(pars.s), equals(ll))
expect_that(lik.s.d(pars.s), equals(ll))
expect_that(lik.s.g(pars.s), equals(ll))
expect_that(lik.s.G(pars.s), equals(ll))

set.seed(1)
pars.s2 <- pars.s + runif(length(pars.s), 0, .1)
ll <- -114.292617388241
expect_that(lik.s.0(pars.s2), equals(ll))
expect_that(lik.s.d(pars.s2), equals(ll, tolerance=1e-7))
expect_that(lik.s.g(pars.s2), equals(ll))
expect_that(lik.s.G(pars.s2), equals(ll))
