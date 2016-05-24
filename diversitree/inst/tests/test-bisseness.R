source("helper-diversitree.R")

## Will be useful to have a function for testing tolerance to within
## 1e-7, as that works out to be how accurate most things actually
## are.
equals7 <- function(...)
  equals(..., tolerance=1e-7)

context("BiSSE-ness")

## First we simulat a 50 species tree, assuming cladogenetic shifts in 
## the trait (i.e., the trait only changes at speciation).
## Red is state '1', black is state '0', and we let red lineages
## speciate at twice the rate of black lineages.
## The simulation starts in state 0.
set.seed(3)
pars <- c(0.1, 0.2, 0.03, 0.03, 0, 0, 0.1, 0, 0.1, 0)
phy <- tree.bisseness(pars, max.taxa=50, x0=0)

## Different control parameters
control.d <- list(backend="deSolve")
control.g <- list(backend="gslode", compiled=FALSE)
control.G <- list(backend="gslode", compiled=TRUE)
control.i <- list(backend="invalid_backend")

lik.0 <- make.bisseness(phy, phy$tip.state)
lik.d <- make.bisseness(phy, phy$tip.state, control=control.d)
lik.g <- make.bisseness(phy, phy$tip.state, control=control.g)
lik.G <- make.bisseness(phy, phy$tip.state, control=control.G)

ll1 <- -174.795395221665
expect_that(lik.0(pars), equals(ll1))
expect_that(lik.d(pars), equals7(ll1))
expect_that(lik.g(pars), equals(ll1))
expect_that(lik.G(pars), equals(ll1))

set.seed(1)
pars2 <- pars + runif(10, 0, .1)

opts <- data.frame(surv=rep(c(TRUE, FALSE), each=3),
                   root=rep(c(ROOT.FLAT, ROOT.OBS, ROOT.EQUI), 2),
                   ll=c(-178.272212778053,
                     -178.264430041395,
                     -178.279034591254,
                     -181.723543144728, 
                     -181.577594926641,
                     -181.836785330049))

liks <- list(lik.0, lik.d, lik.g, lik.G)

## To generate the list above:
## dput(mapply(function(r, s) lik.0(pars2, root=r, condition.surv=s),
##             opts$root, opts$surv))
for ( f in liks )
  for ( i in seq_len(nrow(opts)) )
    expect_that(f(pars2, root=opts$root[i], condition.surv=opts$surv[i]),
                equals7(opts$ll[i]))

## Unresolved tip clade: Here we collapse one clade in the 50 species
## tree (involving sister species sp70 and sp71) and illustrate the use
## of BiSSEness with unresolved tip clades.
phy.u <- drop.tip(phy, "sp71")
states <- phy.u$tip.state[phy.u$tip.label]
states["sp70"] <- NA
unresolved <- data.frame(tip.label=c("sp70"), Nc=2, n0=2, n1=0)

## This builds the likelihood of the data according to BiSSEness:
## SW: warning suppression here because of the warning about
## unresolved clades not being tested extensively.
lik.unresolved <-
  suppressWarnings(make.bisseness(phy.u, states, unresolved))
## e.g., the likelihood of the true parameters is:
expect_that(lik.unresolved(pars), equals(-174.657490373083))
expect_that(lik.unresolved(pars2), equals(-178.182216187058))
