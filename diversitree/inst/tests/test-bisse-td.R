source("helper-diversitree.R")

context("BiSSE (time chunks)")

set.seed(4)
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
phy <- tree.bisse(pars, max.t=30, x0=0)

lik.b <- make.bisse(phy, phy$tip.state)
lik.td <- make.bisse.td(phy, phy$tip.state, 2)

## First, evaluate the functions with no time effect and check that they
## are the same as the base BiSSE model
t <- max(branching.times(phy))/2
p.td <- c(t, pars, pars)

ll <- lik.b(pars)
expect_that(ll, equals(-159.709959436709))
expect_that(lik.td(p.td), equals(ll))

set.seed(1)
pars2 <- pars + runif(6, 0, .2)
p2.td <- c(t, pars2, pars2)  

ll <- -189.608098958571
expect_that(lik.b(pars2), equals(ll))
expect_that(lik.td(p2.td), equals(ll))

set.seed(1)
p3.td <- p2.td + runif(length(p2.td), 0, .2)
expect_that(lik.td(p3.td), equals(-183.125434819243))

lik.t <- make.bisse.t(phy, phy$tip.state, rep("stepf.t", 6))
t <- p3.td[1]
p <- matrix(p3.td[-1], ncol=2)
p3.t <- c(p[1,], t,
          p[2,], t,
          p[3,], t,
          p[4,], t,
          p[5,], t,
          p[6,], t)
names(p3.t) <- argnames(lik.t)
names(p3.td) <- argnames(lik.td)
expect_that(lik.t(p3.t), equals(lik.td(p3.td)))

context("MuSSE (time chunks)")

pars <- c(.1,  .15,  .2,  # lambda 1, 2, 3
          .03, .045, .06, # mu 1, 2, 3
          .05, 0,         # q12, q13
          .05, .05,       # q21, q23
          0,   .05)       # q31, q32
set.seed(2)
phy <- tree.musse(pars, 30, x0=1)

## For comparison, make a plain MuSSE likelihood function
lik.m <- make.musse(phy, phy$tip.state, 3)

## Create the time-dependent likelihood function.  The final argument
## here is the number of 'epochs' that are allowed.  Two epochs is one
## switch point.
lik.td <- make.musse.td(phy, phy$tip.state, 3, 2)

t <- max(branching.times(phy))/2  
p.td <- c(t, pars, pars)
names(p.td) <- argnames(lik.td)

expect_that(lik.m(pars), equals(-110.836441114592))
expect_that(lik.td(p.td), equals(lik.m(pars)))

set.seed(1)
pars2 <- pars + runif(length(pars), 0, .2)
p2.td <- c(t, pars2, pars2)

expect_that(lik.m(pars2), equals(-118.952881268597))
expect_that(lik.td(p2.td), equals(lik.m(pars2)))

set.seed(1)
p3.td <- p2.td + runif(length(p2.td), 0, .2)
expect_that(lik.td(p3.td), equals(-122.524082065168))
