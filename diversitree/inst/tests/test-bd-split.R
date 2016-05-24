## TODO: Could do with more extensive tests here.
source("helper-diversitree.R")

## Will be useful to have a function for testing tolerance to within
## 1e-7, as that works out to be how accurate most things actually
## are.
equals7 <- function(...)
  equals(..., tolerance=1e-7)

context("Birth-death (split)")

set.seed(1)
pars <- c(.1, .03)
phy <- trees(pars, "bd", max.taxa=30)[[1]]

## Construct the plain likelihood function as a benchmark:
lik <- make.bd(phy)
expect_that(lik(pars), equals(-21.7455391151786))

## Split this phylogeny at three points: nd11, nd13 and nd26
nodes <- c("nd11", "nd13", "nd26")

## This is the index in ape's node indexing:
nodes.i <- match(nodes, phy$node.label) + length(phy$tip.label)

## To make a split likelihood function, pass the node locations and times in:
lik.s <- make.bd.split(phy, nodes)

## The parameters must be a list of the same length as the number of
## partitions.  Partition '1' is the root partition, and partition i is
## the partition rooted at the node[i-1]
pars4 <- rep(pars, 4)
names(pars4) <- argnames(lik.s)

## Run the likelihod calculation:
expect_that(lik.s(pars4), equals(lik(pars)))

set.seed(1)
pars2 <- pars + runif(2, 0, .1)
expect_that(lik(pars2), equals(-21.9046823306346))
expect_that(lik.s(rep(pars2, 4)), equals(lik(pars2)))

set.seed(1)
pars4.2 <- pars4 + runif(length(pars4), 0, .1)
expect_that(lik.s(pars4.2), equals(-23.48722743120750))

## You can use the labelled nodes rather than indices:
lik.s2 <- make.bd.split(phy, nodes)
expect_that(lik.s(pars4), is_identical_to(lik.s2(pars4)))

## All the usual ML/MCMC functions work as before:
fit <- find.mle(lik.s, pars4)
expect_that(fit$lnLik, equals(-19.65871005552911))

