source("helper-diversitree.R")

## Will be useful to have a function for testing tolerance to within
## 1e-7, as that works out to be how accurate most things actually
## are.
equals7 <- function(...)
  equals(..., tolerance=1e-7)
equals6 <- function(...)
  equals(..., tolerance=1e-6)

context("GeoSSE")

tree <- read.tree(text="((((0:0.461876,(1:0.307717,(2:0.231825,3:0.231825):0.075892):0.154159):0.425922,((4:0.244819,5:0.244819):0.004749,6:0.249568):0.638231):0.142051,7:1.029850):0.038423,(((8:0.510933,(9:0.427929,(10:0.119778,11:0.119778):0.308151):0.083004):0.007428,(12:0.488316,13:0.488316):0.030044):0.100160,14:0.618521):0.449752);")
states <- c(0, 1, 0, 2, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0)
names(states) <- as.character(seq(Ntip(tree))-1)

rate.names <- c("sA", "sB", "sAB", "xA", "xB", "dA", "dB")
pars0 <- c(0.9, 0.8, 0.1, 0.2, 0.3, 0.5, 0.6)

## Test the starting parameters:
ans <- starting.point.geosse(tree)
ans0 <- c( rep(3.339799, 3), rep(1.669899, 4) )
names(ans0) <- rate.names
expect_that(ans, equals6(ans0))

ans <- starting.point.geosse(tree, eps=0)
ans0 <- c( rep(1.8861312, 3), rep(0, 2), rep(0.1886131, 2) )
names(ans0) <- rate.names
expect_that(ans, equals(ans0))

## Likelihood calculations:

## Different control parameters
control.d <- list(backend="deSolve")
control.g <- list(backend="gslode", compiled=FALSE)
control.G <- list(backend="gslode", compiled=TRUE)
control.i <- list(backend="invalid_backend")


pars <- pars0

lik.0 <- make.geosse(tree, states)
lik.d <- make.geosse(tree, states, control=control.d)
lik.g <- make.geosse(tree, states, control=control.g)
lik.G <- make.geosse(tree, states, control=control.G)

ll <- lik.0(pars)
expect_that(ll, equals(-24.0712795540242))
expect_that(lik.d(pars), equals(ll))
expect_that(lik.g(pars), equals(ll))
expect_that(lik.G(pars), equals(ll))

ll <- lik.0(pars, condition.surv=FALSE)
expect_that(ll, equals(-23.7157388738112))

ll <- lik.0(pars, root=ROOT.EQUI)
expect_that(ll, equals(-24.2550850966392))

ll <- lik.0(pars, root=ROOT.GIVEN, root.p=c(0.6,0.4,0.2))
expect_that(ll, equals(-24.1432006787187))

ll.cmp <- lik.0(pars, root=ROOT.GIVEN, root.p=rep(1,3)/3)
ll <- lik.0(pars, root=ROOT.FLAT)
expect_that(ll, equals(ll.cmp))

lnL <- constrain(lik.0, sAB ~ 0)
ll <- lnL(pars[-3])
expect_that(ll, equals(-23.868263868013))

lnL <- constrain(lik.0, dB ~ dA, sAB ~ 0)
ll <- lnL(pars[-c(3, 7)])
expect_that(ll, equals(-24.0579179253941))

## MLE:
pars <- pars0
names(pars) <- argnames(lik.0)
ans <- find.mle(lik.0, pars)
ans0.par <- c(1.48314580114945, 0.365350237345147, 
              3.17741161756545e-06, 3.4178455229063e-08,
              1.31021569005515e-05, 1.27520486435171,
              1.20617604935305)
names(ans0.par) <- rate.names
expect_that(coef(ans), equals7(ans0.par))
expect_that(logLik(ans)[1], equals(-19.3289142828181))

lik.c <- constrain(lik.0, dB ~ dA, sAB ~ 0)
pars <- pars0[-c(3,7)]
names(pars) <- argnames(lik.c)
ans <- find.mle(lik.c, pars)
ans0.par <- c(1.48147003955874, 0.369816921968662, 
              9.21189636011179e-09, 1.52697977801869e-08,
              1.27032432170033)
names(ans0.par) <- rate.names[-c(3,7)]
expect_that(coef(ans), equals(ans0.par))
expect_that(logLik(ans)[1], equals(-19.3290833415077))

## Old MCMC test:

## set.seed(1)
## pars <- pars0
## ans <- mcmc(lnL.full, pars, nsteps=3, lower=0, upper=3, w=3/5, 
##             prior=make.prior.exponential(seq(0.8, by=0.1, length.out=7)),
##             print.every=0) 
## checkEquals(ans$p[3], -30.56511, tolerance=tol)

## set.seed(1)
## lnL <- constrain(lnL.full, dB ~ dA, sAB ~ 0)
## pars <- pars0[-c(3,7)]
## ans <- mcmc(lnL, pars, nsteps=3, lower=0, upper=3, w=seq(0.8, by=0.1,
##             length.out=5), prior=make.prior.exponential(1), print.every=0)
## checkEquals(ans$p[3], -28.18477, tolerance=tol)


## Old Split test:
## with 0.8-4, only split.t = Inf works (for any model)
## lnL.split <- make.geosse.split(tree, states, c(27, 29), split.t = Inf)
## checkEquals(lnL.full(pars0), lnL.split(rep(pars0, 3)), tolerance=tol)
## checkEquals(lnL.split(c(pars0, pars0*1.5, pars0*0.5)), -23.82277, tolerance=tol)

pars.g <- pars0
names(pars.g) <- diversitree:::default.argnames.geosse()
pars.c <- diversitree:::pars.ge.to.cl(pars.g)

set.seed(1)
phy <- trees(pars.g, type="geosse", n=2, max.t=4, x0=0)
expect_that(lapply(phy, Ntip), equals(list(20, 132)))
lnL <- make.geosse(phy[[2]], phy[[2]]$tip.state)
expect_that(lnL(pars.g), equals(-252.717350358088))

## Does this belong in ClaSSE's tests?
set.seed(3)
phy2 <- tree.classe(pars.c, max.t=5)
lnL.c <- make.classe(phy2, phy2$tip.state, 3)
lnL.g <- make.geosse(phy2, phy2$tip.state-1)
expect_that(lnL.c(pars.c), equals(lnL.g(pars.g)))
