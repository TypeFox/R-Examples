source("helper-diversitree.R")

context("QuaSSE (split)")

## Build tree:
lambda <- function(x) sigmoid.x(x, 0.1, 0.2,  0, 2.5)
mu <- function(x) constant.x(x, 0.03)
char <- make.brownian.with.drift(0, 0.025)

## set.seed(1)
## phy <- tree.quasse(c(lambda, mu, char), max.taxa=15, x0=0,
##                    single.lineage=FALSE, verbose=TRUE)
load("phy.Rdata")

pars <- c(.1, .2, 0, 2.5, .03, 0, .01)
pars.s <- rep(pars, 2)  

sd <- 1/200
control.C.1 <- list(dt.max=1/200)
control.C.2 <- c(control.C.1, tips.combined=TRUE)
control.M.1 <- list(method="mol")
control.R.1 <- list(dt.max=1/200, method="fftR")

if (check.fftC(FALSE)) {
lik.s <- make.quasse.split(phy, phy$tip.state, sd, sigmoid.x,
                           constant.x, "nd5", Inf, control.C.1)
lik.q <- make.quasse(phy, phy$tip.state, sd, sigmoid.x, constant.x,
                     control.C.1)
ll.q <- lik.q(pars)
expect_that(ll.q, equals(-62.06409424693976))

pars.s <- rep(pars, 2)
names(pars.s) <- argnames(lik.s)
expect_that(lik.s(pars.s), equals(ll.q))

set.seed(1)
pars2 <- pars + runif(length(pars), 0, .05)
pars2.s <- rep(pars2, 2)
ll.q <- lik.q(pars2)
expect_that(ll.q, equals(-55.67237675384200))
expect_that(lik.s(pars2.s), equals(ll.q))

pars3.s <- pars + runif(length(pars.s), 0, .05)
expect_that(lik.s(pars3.s), equals(-54.47383577050427))
}
