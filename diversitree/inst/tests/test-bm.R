source("helper-diversitree.R")

no.stdout <- function(expr) {
  sink(tempfile())
  on.exit(sink())
  eval.parent(substitute(expr))
}

context("Brownian motion")

## Simulated tree and traits:
set.seed(1)
phy <- tree.bd(pars=c(1,0), max.taxa=100)
states <- sim.character(phy, 1)
se <- 0.1

## Currently three different calculators, and two different backends.
## We may get an alternative backend for contrasts soon.
lik.vcv <- make.bm(phy, states, control=list(method="vcv"))
lik.pru.R <- make.bm(phy, states,
                     control=list(method="pruning", backend="R"))
lik.pru.C <- make.bm(phy, states,
                     control=list(method="pruning", backend="C"))
lik.con <- make.bm(phy, states, control=list(method="contrasts"))

lik.vcv.se <- make.bm(phy, states, se, control=list(method="vcv"))
lik.pru.R.se <- make.bm(phy, states, se,
                     control=list(method="pruning", backend="R"))
lik.pru.C.se <- make.bm(phy, states, se,
                     control=list(method="pruning", backend="C"))
## Not yet supported
expect_that(make.bm(phy, states, se,
                    control=list(method="contrasts")),
            throws_error())

## First, a simple test on a value that is known.  Obviously if the
## tree simulators change this will break!  But all four will then
## also break.
test_that("Likelihood calculations agree on known case", {
  ll <- -128.150053529354
  expect_that(lik.vcv(1), equals(ll))
  expect_that(lik.pru.R(1), equals(ll))
  expect_that(lik.pru.C(1), equals(ll))
  expect_that(lik.con(1), equals(ll))
})

test_that("Likelihood calculations agree on known case", {
  ll <- -128.150053529354
  expect_that(lik.vcv(1), equals(ll))
  expect_that(lik.pru.R(1), equals(ll))
  expect_that(lik.pru.C(1), equals(ll))
  expect_that(lik.con(1), equals(ll))
})


## Exercise the different root treatents.
s2 <- 0.9 # diffusion parameter for all cases below

test_that("Calculations agree for ROOT.MAX", {
  ll.max  <- -127.97762790984
  expect_that(lik.vcv(s2, root=ROOT.MAX), equals(ll.max))
  expect_that(lik.pru.R(s2, root=ROOT.MAX), equals(ll.max))
  expect_that(lik.pru.C(s2, root=ROOT.MAX), equals(ll.max))
  expect_that(lik.con(s2, root=ROOT.MAX), equals(ll.max))
})

test_that("Calculations agree for ROOT.FLAT", {
  ll.flat <- -127.507880373793
  expect_that(lik.vcv(s2, root=ROOT.FLAT), throws_error())
  expect_that(lik.pru.R(s2, root=ROOT.FLAT), equals(ll.flat))
  expect_that(lik.pru.C(s2, root=ROOT.FLAT), equals(ll.flat))
  expect_that(lik.con(s2, root=ROOT.FLAT), equals(ll.flat))
})

test_that("Calculations agree for ROOT.OBS", {
  ll.obs  <- -128.32420150012
  expect_that(lik.vcv(s2, root=ROOT.OBS), throws_error())
  expect_that(lik.pru.R(s2, root=ROOT.OBS), equals(ll.obs))
  expect_that(lik.pru.C(s2, root=ROOT.OBS), equals(ll.obs))
  expect_that(lik.con(s2, root=ROOT.OBS), equals(ll.obs))
})

test_that("Calculations agree for ROOT.GIVEN", {
  ll.0    <- -127.999550052743 # root.x = 0
  expect_that(lik.vcv(s2, root=ROOT.GIVEN, root.x=0), throws_error())
  expect_that(lik.pru.R(s2, root=ROOT.GIVEN, root.x=0), equals(ll.0))
  expect_that(lik.pru.C(s2, root=ROOT.GIVEN, root.x=0), equals(ll.0))
  expect_that(lik.con(s2, root=ROOT.GIVEN, root.x=0), equals(ll.0))

  ll.1    <- -129.55548714521  # root.x = 1
  expect_that(lik.vcv(s2, root=ROOT.GIVEN, root.x=1), throws_error())
  expect_that(lik.pru.R(s2, root=ROOT.GIVEN, root.x=1), equals(ll.1))
  expect_that(lik.pru.C(s2, root=ROOT.GIVEN, root.x=1), equals(ll.1))
  expect_that(lik.con(s2, root=ROOT.GIVEN, root.x=1), equals(ll.1))
})

## Repeat with SE
test_that("Calculations agree for ROOT.MAX (with SE)", {
  ll.max  <- -128.536896031852
  expect_that(lik.vcv.se(s2, root=ROOT.MAX), equals(ll.max))
  expect_that(lik.pru.R.se(s2, root=ROOT.MAX), equals(ll.max))
  expect_that(lik.pru.C.se(s2, root=ROOT.MAX), equals(ll.max))
  expect_that(lik.con.se(s2, root=ROOT.MAX), throws_error())
})

test_that("Calculations agree for ROOT.FLAT (with SE)", {
  ll.flat <- -128.066836613671
  expect_that(lik.vcv.se(s2, root=ROOT.FLAT), throws_error())
  expect_that(lik.pru.R.se(s2, root=ROOT.FLAT), equals(ll.flat))
  expect_that(lik.pru.C.se(s2, root=ROOT.FLAT), equals(ll.flat))
  expect_that(lik.con.se(s2, root=ROOT.FLAT), throws_error())
})

test_that("Calculations agree for ROOT.OBS (with SE)", {
  ll.obs  <- -128.883469622132
  expect_that(lik.vcv.se(s2, root=ROOT.OBS), throws_error())
  expect_that(lik.pru.R.se(s2, root=ROOT.OBS), equals(ll.obs))
  expect_that(lik.pru.C.se(s2, root=ROOT.OBS), equals(ll.obs))
  expect_that(lik.con.se(s2, root=ROOT.OBS), throws_error())
})

test_that("Calculations agree for ROOT.GIVEN (with SE)", {
  ll.0    <- -128.558797137454 # root.x = 0
  expect_that(lik.vcv.se(s2, root=ROOT.GIVEN, root.x=0), throws_error())
  expect_that(lik.pru.R.se(s2, root=ROOT.GIVEN, root.x=0), equals(ll.0))
  expect_that(lik.pru.C.se(s2, root=ROOT.GIVEN, root.x=0), equals(ll.0))
  expect_that(lik.con.se(s2, root=ROOT.GIVEN, root.x=0), throws_error())

  ll.1    <- -130.113708854126  # root.x = 1
  expect_that(lik.vcv.se(s2, root=ROOT.GIVEN, root.x=1), throws_error())
  expect_that(lik.pru.R.se(s2, root=ROOT.GIVEN, root.x=1), equals(ll.1))
  expect_that(lik.pru.C.se(s2, root=ROOT.GIVEN, root.x=1), equals(ll.1))
  expect_that(lik.con.se(s2, root=ROOT.GIVEN, root.x=1), throws_error())
})

## Now, fit models using ML:
test_that("Can fit models with ML", {
  fit.vcv   <- find.mle(lik.vcv,   .1)
  fit.pru.R <- find.mle(lik.pru.R, .1)
  fit.pru.C <- find.mle(lik.pru.C, .1)
  fit.con   <- find.mle(lik.con,   .1)

  expect_that(fit.pru.R, equals(fit.vcv))
  expect_that(fit.pru.C, equals(fit.vcv))
  expect_that(fit.con,   equals(fit.vcv))

  x <- structure(as.numeric(states), names=names(states))
  fit.geiger <- no.stdout(fitContinuous(phy, x,
                                        control=list(niter=5)))
  expect_that(fit.vcv$lnLik, equals(fit.geiger$opt$lnL))
  expect_that(fit.vcv$par, equals(fit.geiger$opt$sigsq,
                                  tolerance=1e-4,
                                  check.attributes=FALSE))
})

test_that("Can fit models with ML (with SE)", {
  fit.vcv   <- find.mle(lik.vcv.se,   .1)
  fit.pru.R <- find.mle(lik.pru.R.se, .1)
  fit.pru.C <- find.mle(lik.pru.C.se, .1)
  # fit.con   <- find.mle(lik.con,   .1)

  expect_that(fit.pru.R, equals(fit.vcv))
  expect_that(fit.pru.C, equals(fit.vcv))
  # expect_that(fit.con,   equals(fit.vcv))

  x <- structure(as.numeric(states), names=names(states))
  fit.geiger <- no.stdout(fitContinuous(phy, x, SE=se,
                                        control=list(niter=5)))
  expect_that(fit.vcv$lnLik, equals(fit.geiger$opt$lnL))
  expect_that(fit.vcv$par, equals(fit.geiger$opt$sigsq,
                                  tolerance=1e-4,
                                  check.attributes=FALSE))
})
