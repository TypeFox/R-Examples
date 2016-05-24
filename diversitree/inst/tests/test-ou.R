source("helper-diversitree.R")

no.stdout <- function(expr) {
  sink(tempfile())
  on.exit(sink())
  eval.parent(substitute(expr))
}

context("Ornstein-Uhlenbeck")

## Simulated tree and traits:
set.seed(1)
phy <- tree.bd(pars=c(1,0), max.taxa=100)
states <- sim.character(phy, 1)
se <- 0.1

## First, without optimum (this works for all methods).
lik.vcv <- make.ou(phy, states, control=list(method="vcv"))
lik.pru.R <- make.ou(phy, states,
                     control=list(method="pruning", backend="R"))
lik.pru.C <- make.ou(phy, states,
                     control=list(method="pruning", backend="C"))
lik.con <- make.ou(phy, states, control=list(method="contrasts"))

lik.vcv.se <- make.ou(phy, states, se, control=list(method="vcv"))
lik.pru.R.se <- make.ou(phy, states, se,
                     control=list(method="pruning", backend="R"))
lik.pru.C.se <- make.ou(phy, states, se,
                     control=list(method="pruning", backend="C"))
## Not yet supported
expect_that(make.ou(phy, states, se,
                    control=list(method="contrasts")),
            throws_error())

## First, a simple test on a value that is known.  Obviously if the
## tree simulators change this will break!  But all four will then
## also break.
test_that("Likelihood calculations agree on known case", {
  ## Start with BM; set alpha to zero:
  pars <- c(1, 0)
  ll <- -128.150053529354
  expect_that(lik.vcv(pars), equals(ll))
  expect_that(lik.pru.R(pars), equals(ll))
  expect_that(lik.pru.C(pars), equals(ll))
  expect_that(lik.con(pars), equals(ll))

  ## Next, bump up alpha a little.
  pars <- c(1, .1)
  ll <- -127.649717292364
  expect_that(lik.vcv(pars), equals(ll))
  expect_that(lik.pru.R(pars), equals(ll))
  expect_that(lik.pru.C(pars), equals(ll))
  expect_that(lik.con(pars), equals(ll))

  ## Calculations can become unstable with |alpha| << 1, so test that
  ## we get the BM likelihood there with extremely small alpha values,
  ## rathe than just bailing.
  pars <- c(1, 1e-20)
  ll <- -128.150053529354
  expect_that(lik.vcv(pars), equals(ll))
  expect_that(lik.pru.R(pars), equals(ll))
  expect_that(lik.pru.C(pars), equals(ll))
  expect_that(lik.con(pars), equals(ll))
})

## Exercise the different root treatents.
pars <- c(1, .1) # parameters for all cases below

test_that("Calculations agree for ROOT.MAX", {
  ll.max  <- -127.649717292364
  expect_that(lik.vcv(pars, root=ROOT.MAX), equals(ll.max))
  expect_that(lik.pru.R(pars, root=ROOT.MAX), equals(ll.max))
  expect_that(lik.pru.C(pars, root=ROOT.MAX), equals(ll.max))
  expect_that(lik.con(pars, root=ROOT.MAX), equals(ll.max))
})

test_that("Calculations agree for ROOT.FLAT", {
  ll.flat  <- -127.382402189285
  expect_that(lik.vcv(pars, root=ROOT.FLAT), throws_error())
  expect_that(lik.pru.R(pars, root=ROOT.FLAT), equals(ll.flat))
  expect_that(lik.pru.C(pars, root=ROOT.FLAT), equals(ll.flat))
  expect_that(lik.con(pars, root=ROOT.FLAT), equals(ll.flat))
})

test_that("Calculations agree for ROOT.OBS", {
  ll.obs  <- -127.996290882644
  expect_that(lik.vcv(pars, root=ROOT.OBS), throws_error())
  expect_that(lik.pru.R(pars, root=ROOT.OBS), equals(ll.obs))
  expect_that(lik.pru.C(pars, root=ROOT.OBS), equals(ll.obs))
  expect_that(lik.con(pars, root=ROOT.OBS), equals(ll.obs))
})

test_that("Calculations agree for ROOT.GIVEN", {
  ll.0    <- -127.679069119391 # root.x = 0
  expect_that(lik.vcv(pars, root=ROOT.GIVEN, root.x=0), throws_error())
  expect_that(lik.pru.R(pars, root=ROOT.GIVEN, root.x=0), equals(ll.0))
  expect_that(lik.pru.C(pars, root=ROOT.GIVEN, root.x=0), equals(ll.0))
  expect_that(lik.con(pars, root=ROOT.GIVEN, root.x=0), equals(ll.0))

  ll.1    <- -129.984551312729  # root.x = 1
  expect_that(lik.vcv(pars, root=ROOT.GIVEN, root.x=1), throws_error())
  expect_that(lik.pru.R(pars, root=ROOT.GIVEN, root.x=1), equals(ll.1))
  expect_that(lik.pru.C(pars, root=ROOT.GIVEN, root.x=1), equals(ll.1))
  expect_that(lik.con(pars, root=ROOT.GIVEN, root.x=1), equals(ll.1))
})

## Again, with std error
test_that("Calculations agree for ROOT.MAX (with SE)", {
  ll.max  <- -128.296638672837
  expect_that(lik.vcv.se(pars, root=ROOT.MAX), equals(ll.max))
  expect_that(lik.pru.R.se(pars, root=ROOT.MAX), equals(ll.max))
  expect_that(lik.pru.C.se(pars, root=ROOT.MAX), equals(ll.max))
  expect_that(lik.con.se(pars, root=ROOT.MAX), throws_error())
})

test_that("Calculations agree for ROOT.FLAT (with SE)", {
  ll.flat <- -128.028893988585
  expect_that(lik.vcv.se(pars, root=ROOT.FLAT), throws_error())
  expect_that(lik.pru.R.se(pars, root=ROOT.FLAT), equals(ll.flat))
  expect_that(lik.pru.C.se(pars, root=ROOT.FLAT), equals(ll.flat))
  expect_that(lik.con.se(pars, root=ROOT.FLAT), throws_error())
})

test_that("Calculations agree for ROOT.OBS (with SE)", {
  ll.obs  <- -128.643212263117
  expect_that(lik.vcv.se(pars, root=ROOT.OBS), throws_error())
  expect_that(lik.pru.R.se(pars, root=ROOT.OBS), equals(ll.obs))
  expect_that(lik.pru.C.se(pars, root=ROOT.OBS), equals(ll.obs))
  expect_that(lik.con.se(pars, root=ROOT.OBS), throws_error())
})

test_that("Calculations agree for ROOT.GIVEN (with SE)", {
  ll.0    <- -128.325962379357
  expect_that(lik.vcv.se(pars, root=ROOT.GIVEN, root.x=0), throws_error())
  expect_that(lik.pru.R.se(pars, root=ROOT.GIVEN, root.x=0), equals(ll.0))
  expect_that(lik.pru.C.se(pars, root=ROOT.GIVEN, root.x=0), equals(ll.0))
  expect_that(lik.con.se(pars, root=ROOT.GIVEN, root.x=0), throws_error())

  ll.1    <- -130.629441568798
  expect_that(lik.vcv.se(pars, root=ROOT.GIVEN, root.x=1), throws_error())
  expect_that(lik.pru.R.se(pars, root=ROOT.GIVEN, root.x=1), equals(ll.1))
  expect_that(lik.pru.C.se(pars, root=ROOT.GIVEN, root.x=1), equals(ll.1))
  expect_that(lik.con.se(pars, root=ROOT.GIVEN, root.x=1), throws_error())
})

## Now, fit models using ML:
test_that("Can fit models with ML", {
  fit.vcv   <- find.mle(lik.vcv,   pars)
  fit.pru.R <- find.mle(lik.pru.R, pars)
  fit.pru.C <- find.mle(lik.pru.C, pars)
  fit.con   <- find.mle(lik.con,   pars)

  expect_that(fit.pru.R, equals(fit.vcv))
  expect_that(fit.pru.C, equals(fit.vcv))
  expect_that(fit.con,   equals(fit.vcv))

  x <- structure(as.numeric(states), names=names(states))
  fit.geiger <- no.stdout(fitContinuous(phy, x, model="OU",
                                        control=list(niter=5)))
  expect_that(fit.vcv$lnLik, equals(fit.geiger$opt$lnL))
  expect_that(coef(fit.vcv),
              equals(unlist(fit.geiger$opt[c("sigsq", "alpha")]),
                                  tolerance=1e-4,
                                  check.attributes=FALSE))
})

test_that("Can fit models with ML (with SE)", {
  fit.vcv   <- find.mle(lik.vcv.se,   pars)
  fit.pru.R <- find.mle(lik.pru.R.se, pars)
  fit.pru.C <- find.mle(lik.pru.C.se, pars)
  # fit.con   <- find.mle(lik.con.se,   pars)

  expect_that(fit.pru.R, equals(fit.vcv))
  expect_that(fit.pru.C, equals(fit.vcv))
  # expect_that(fit.con,   equals(fit.vcv))

  x <- structure(as.numeric(states), names=names(states))
  fit.geiger <- no.stdout(fitContinuous(phy, x, SE=se, model="OU",
                                        control=list(niter=5)))
  expect_that(fit.vcv$lnLik, equals(fit.geiger$opt$lnL))
  expect_that(coef(fit.vcv),
              equals(unlist(fit.geiger$opt[c("sigsq", "alpha")]),
                                  tolerance=1e-4,
                                  check.attributes=FALSE))
})

## Now, continue on with the "with-optimum" version:
expect_that(lik.vcv.opt <- make.ou(phy, states, with.optimum=TRUE,
                                   control=list(method="vcv")),
            throws_error())
lik.pru.R.opt <- make.ou(phy, states, with.optimum=TRUE,
                         control=list(method="pruning", backend="R"))
lik.pru.C.opt <- make.ou(phy, states, with.optimum=TRUE,
                         control=list(method="pruning", backend="C"))
expect_that(lik.con.opt <- make.ou(phy, states, with.optimum=TRUE,
                                   control=list(method="contrasts")),
            throws_error())

expect_that(lik.vcv.se.opt <- make.ou(phy, states, se,
                                      with.optimum=TRUE,
                                      control=list(method="vcv")),
            throws_error())
lik.pru.R.se.opt <- make.ou(phy, states, se, with.optimum=TRUE,
                        control=list(method="pruning", backend="R"))
lik.pru.C.se.opt <- make.ou(phy, states, se, with.optimum=TRUE,
                     control=list(method="pruning", backend="C"))
expect_that(make.ou(phy, states, se, with.optimum=TRUE,
                    control=list(method="contrasts")),
            throws_error())

## at the moment, the tests are pretty straightforward and not
## comprehensive.  The differences are not very great between the
## with-optimum and without-optimum cases are actually fairly small,
## so nailing these down is tricky.
test_that("Likelihood calculations agree on known case", {
  ## Start with BM; set alpha to zero:
  pars <- c(1, 0, 0)
  ll <- -128.150053529354
  expect_that(lik.pru.R.opt(pars), equals(ll))
  expect_that(lik.pru.C.opt(pars), equals(ll))

  ## Next, bump up alpha a little.
  pars <- c(1, .1, 0)
  ll <- -127.649717292364
  expect_that(lik.pru.R.opt(pars), equals(ll, tolerance=1e-10))
  expect_that(lik.pru.C.opt(pars), equals(ll, tolerance=1e-10))

  ## In this parameter space, we barely depend on theta, but should
  ## see some change with absolutely massive changes in theta:
  pars2 <- pars + c(0, 0, 10000)
  expect_that(abs(lik.pru.R.opt(pars) - lik.pru.R.opt(pars2)) > 0,
              is_true())
  expect_that(abs(lik.pru.C.opt(pars) - lik.pru.C.opt(pars2)) > 0,
              is_true())

  ## Calculations can become unstable with |alpha| << 1, so test that
  ## we get the BM likelihood there with extremely small alpha values,
  ## rathe than just bailing.
  pars <- c(1, 1e-20, 0)
  ll <- -128.150053529354
  expect_that(lik.pru.R.opt(pars), equals(ll))
  expect_that(lik.pru.C.opt(pars), equals(ll))
})

pars <- c(1, .1, 0)
test_that("Calculations agree for ROOT.MAX", {
  ll.max  <- -127.649717292364
  expect_that(lik.pru.R.opt(pars, root=ROOT.MAX), equals(ll.max))
  expect_that(lik.pru.C.opt(pars, root=ROOT.MAX), equals(ll.max))
})

test_that("Calculations agree for ROOT.FLAT", {
  ll.flat <- -127.382402189285 # no opt
  ll.flat <- -127.012706241672
  expect_that(lik.pru.R.opt(pars, root=ROOT.FLAT), equals(ll.flat))
  expect_that(lik.pru.C.opt(pars, root=ROOT.FLAT), equals(ll.flat))
})

test_that("Calculations agree for ROOT.OBS", {
  ll.obs <- -127.996290882644 # no opt
  ll.obs <- -127.996290882644
  expect_that(lik.pru.R.opt(pars, root=ROOT.OBS), equals(ll.obs))
  expect_that(lik.pru.C.opt(pars, root=ROOT.OBS), equals(ll.obs))
})

test_that("Calculations agree for ROOT.GIVEN", {
  ll.0 <- -127.679069119391 # root.x = 0, no opt
  ll.0 <- -127.679069119391 # root.x = 0
  expect_that(lik.pru.R.opt(pars, root=ROOT.GIVEN, root.x=0), equals(ll.0))
  expect_that(lik.pru.C.opt(pars, root=ROOT.GIVEN, root.x=0), equals(ll.0))

  ll.1 <- -129.984551312729 # root.x = 1, no opt
  ll.1 <- -128.878983743877 # root.x = 1
  expect_that(lik.pru.R.opt(pars, root=ROOT.GIVEN, root.x=1), equals(ll.1))
  expect_that(lik.pru.C.opt(pars, root=ROOT.GIVEN, root.x=1), equals(ll.1))
})

## Again, with std error
test_that("Calculations agree for ROOT.MAX (with SE)", {
  ll.max  <- -128.296638672837
  expect_that(lik.pru.R.se.opt(pars, root=ROOT.MAX), equals(ll.max))
  expect_that(lik.pru.C.se.opt(pars, root=ROOT.MAX), equals(ll.max))
})

test_that("Calculations agree for ROOT.FLAT (with SE)", {
  ll.flat <- -128.028893988585 # no opt
  ll.flat <- -127.659198040972
  expect_that(lik.pru.R.se.opt(pars, root=ROOT.FLAT), equals(ll.flat))
  expect_that(lik.pru.C.se.opt(pars, root=ROOT.FLAT), equals(ll.flat))
})

test_that("Calculations agree for ROOT.OBS (with SE)", {
  ll.obs <- -128.643212263117 # no opt
  ll.obs <- -128.643212263117
  expect_that(lik.pru.R.se.opt(pars, root=ROOT.OBS), equals(ll.obs))
  expect_that(lik.pru.C.se.opt(pars, root=ROOT.OBS), equals(ll.obs))
})

test_that("Calculations agree for ROOT.GIVEN (with SE)", {
  ll.0    <- -128.325962379357 # no opt
  ll.0    <- -128.325962379357
  expect_that(lik.pru.R.se.opt(pars, root=ROOT.GIVEN, root.x=0), equals(ll.0))
  expect_that(lik.pru.C.se.opt(pars, root=ROOT.GIVEN, root.x=0), equals(ll.0))

  ll.1 <- -130.629441568798
  ll.1 <- -129.52483058429
  expect_that(lik.pru.R.se.opt(pars, root=ROOT.GIVEN, root.x=1), equals(ll.1))
  expect_that(lik.pru.C.se.opt(pars, root=ROOT.GIVEN, root.x=1), equals(ll.1))
})
