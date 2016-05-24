source("helper-diversitree.R")

no.stdout <- function(expr) {
  sink(tempfile())
  on.exit(sink())
  eval.parent(substitute(expr))
}

context("Pagel's lambda")

## Simulated tree and traits:
set.seed(1)
phy <- tree.bd(pars=c(1,0), max.taxa=100)
states <- sim.character(phy, 1)
se <- 0.1

## First, without optimum (this works for all methods).
lik.vcv <- make.lambda(phy, states, control=list(method="vcv"))
lik.pru.R <- make.lambda(phy, states,
                     control=list(method="pruning", backend="R"))
lik.pru.C <- make.lambda(phy, states,
                     control=list(method="pruning", backend="C"))
lik.con <- make.lambda(phy, states, control=list(method="contrasts"))

lik.vcv.se <- make.lambda(phy, states, se, control=list(method="vcv"))
lik.pru.R.se <- make.lambda(phy, states, se,
                     control=list(method="pruning", backend="R"))
lik.pru.C.se <- make.lambda(phy, states, se,
                     control=list(method="pruning", backend="C"))
## Not yet supported
expect_that(make.lambda(phy, states, se,
                        control=list(method="contrasts")),
            throws_error())

test_that("Likelihood calculations agree on known case", {
  ## Start with BM; set lambda to one:
  pars <- c(1, 1)
  ll <- -128.150053529354
  expect_that(lik.vcv(pars), equals(ll))
  expect_that(lik.pru.R(pars), equals(ll))
  expect_that(lik.pru.C(pars), equals(ll))
  # expect_that(lik.con(pars), equals(ll))

  ## Next, turn town lambda a little.
  pars <- c(1, 0.7)
  ll <- -156.571120596637

  phy2 <- diversitree:::make.rescale.phylo.lambda(phy)(pars[[2]])
  ll.cmp <- make.bm(phy2, states,
                    control=list(method="vcv"))(pars[[1]])
  expect_that(ll, equals(ll.cmp))
  
  expect_that(lik.vcv(pars), equals(ll))
  expect_that(lik.pru.R(pars), equals(ll))
  expect_that(lik.pru.C(pars), equals(ll))
})

## Excercise the different root treatment.
pars <- c(0.6, 0.7) # parameters for all cases below

test_that("Calculations agree for ROOT.MAX", {
  ll.max  <- -148.105836357897
  expect_that(lik.vcv(pars, root=ROOT.MAX), equals(ll.max))
  expect_that(lik.pru.R(pars, root=ROOT.MAX), equals(ll.max))
  expect_that(lik.pru.C(pars, root=ROOT.MAX), equals(ll.max))
  # expect_that(lik.con(pars, root=ROOT.MAX), equals(ll.max))
})

test_that("Calculations agree for ROOT.FLAT", {
  ll.flat  <- -147.980581724904
  expect_that(lik.vcv(pars, root=ROOT.FLAT), throws_error())
  expect_that(lik.pru.R(pars, root=ROOT.FLAT), equals(ll.flat))
  expect_that(lik.pru.C(pars, root=ROOT.FLAT), equals(ll.flat))
  # expect_that(lik.con(pars, root=ROOT.FLAT), equals(ll.flat))
})

test_that("Calculations agree for ROOT.OBS", {
  ll.obs  <- -148.452409948177
  expect_that(lik.vcv(pars, root=ROOT.OBS), throws_error())
  expect_that(lik.pru.R(pars, root=ROOT.OBS), equals(ll.obs))
  expect_that(lik.pru.C(pars, root=ROOT.OBS), equals(ll.obs))
  # expect_that(lik.con(pars, root=ROOT.OBS), equals(ll.obs))
})

test_that("Calculations agree for ROOT.GIVEN", {
  ll.0    <- -148.14051702295 # root.x = 0
  expect_that(lik.vcv(pars, root=ROOT.GIVEN, root.x=0), throws_error())
  expect_that(lik.pru.R(pars, root=ROOT.GIVEN, root.x=0), equals(ll.0))
  expect_that(lik.pru.C(pars, root=ROOT.GIVEN, root.x=0), equals(ll.0))
  # expect_that(lik.con(pars, root=ROOT.GIVEN, root.x=0), equals(ll.0))

  ll.1    <- -151.168386219595  # root.x = 1
  expect_that(lik.vcv(pars, root=ROOT.GIVEN, root.x=1), throws_error())
  expect_that(lik.pru.R(pars, root=ROOT.GIVEN, root.x=1), equals(ll.1))
  expect_that(lik.pru.C(pars, root=ROOT.GIVEN, root.x=1), equals(ll.1))
  # expect_that(lik.con(pars, root=ROOT.GIVEN, root.x=1), equals(ll.1))
})

## Again, with std error
test_that("Calculations agree for ROOT.MAX (with SE)", {
  ll.max  <- -148.299722874325
  expect_that(lik.vcv.se(pars, root=ROOT.MAX), equals(ll.max))
  expect_that(lik.pru.R.se(pars, root=ROOT.MAX), equals(ll.max))
  expect_that(lik.pru.C.se(pars, root=ROOT.MAX), equals(ll.max))
  # expect_that(lik.con.se(pars, root=ROOT.MAX), throws_error())
})

test_that("Calculations agree for ROOT.FLAT (with SE)", {
  ll.flat <- -148.174001970594
  expect_that(lik.vcv.se(pars, root=ROOT.FLAT), throws_error())
  expect_that(lik.pru.R.se(pars, root=ROOT.FLAT), equals(ll.flat))
  expect_that(lik.pru.C.se(pars, root=ROOT.FLAT), equals(ll.flat))
  # expect_that(lik.con.se(pars, root=ROOT.FLAT), throws_error())
})

test_that("Calculations agree for ROOT.OBS (with SE)", {
  ll.obs  <- -148.646296464604
  expect_that(lik.vcv.se(pars, root=ROOT.OBS), throws_error())
  expect_that(lik.pru.R.se(pars, root=ROOT.OBS), equals(ll.obs))
  expect_that(lik.pru.C.se(pars, root=ROOT.OBS), equals(ll.obs))
  # expect_that(lik.con.se(pars, root=ROOT.OBS), throws_error())
})

test_that("Calculations agree for ROOT.GIVEN (with SE)", {
  ll.0    <- -148.334247766823
  expect_that(lik.vcv.se(pars, root=ROOT.GIVEN, root.x=0), throws_error())
  expect_that(lik.pru.R.se(pars, root=ROOT.GIVEN, root.x=0), equals(ll.0))
  expect_that(lik.pru.C.se(pars, root=ROOT.GIVEN, root.x=0), equals(ll.0))
  # expect_that(lik.con.se(pars, root=ROOT.GIVEN, root.x=0), throws_error())

  ll.1    <- -151.358257137935
  expect_that(lik.vcv.se(pars, root=ROOT.GIVEN, root.x=1), throws_error())
  expect_that(lik.pru.R.se(pars, root=ROOT.GIVEN, root.x=1), equals(ll.1))
  expect_that(lik.pru.C.se(pars, root=ROOT.GIVEN, root.x=1), equals(ll.1))
  # expect_that(lik.con.se(pars, root=ROOT.GIVEN, root.x=1), throws_error())
})

## Now, fit models using ML:
test_that("Can fit models with ML", {
  pars <- c(1, 1)
  fit.vcv   <- find.mle(lik.vcv,   pars)
  fit.pru.R <- find.mle(lik.pru.R, pars)
  fit.pru.C <- find.mle(lik.pru.C, pars)
  ## fit.con   <- find.mle(lik.con,   pars)

  expect_that(fit.pru.R, equals(fit.vcv))
  expect_that(fit.pru.C, equals(fit.vcv))
  ## expect_that(fit.con,   equals(fit.vcv))

  ## Geiger does a surprisingly terrible job of fitting the model,
  ## finding a ML point that is way worse than the global optimum.
  x <- structure(as.numeric(states), names=names(states))
  fit.geiger <-
    no.stdout(suppressWarnings(fitContinuous(phy, x, model="lambda",
                                             control=list(niter=5))))

  ## Even though we don't agree on where the ML point is, let's check
  ## that the likelihood calculations do agree at the respective
  ## points.  First, at the diversitree ML point:
  p <- coef(fit.vcv)[c("lambda", "s2")]
  names(p)[2] <- "sigsq"
  expect_that(as.numeric(fit.geiger$lik(p)),
              equals(fit.vcv$lnLik))

  ## Second, at the geiger ML point.
  p <- coef(fit.geiger)[c("sigsq", "lambda")]
  expect_that(lik.vcv(p), equals(fit.geiger$opt$lnL))
  expect_that(lik.pru.R(p), equals(fit.geiger$opt$lnL))
  expect_that(lik.pru.C(p), equals(fit.geiger$opt$lnL))
})

test_that("Can fit models with ML (with SE)", {
  pars <- c(1, 1)
  fit.vcv   <- find.mle(lik.vcv.se,   pars)
  fit.pru.R <- find.mle(lik.pru.R.se, pars)
  fit.pru.C <- find.mle(lik.pru.C.se, pars)
  ## fit.con   <- find.mle(lik.con.se,   pars)

  expect_that(fit.pru.R, equals(fit.vcv))
  expect_that(fit.pru.C, equals(fit.vcv))
  ## expect_that(fit.con,   equals(fit.vcv))

  x <- structure(as.numeric(states), names=names(states))
  fit.geiger <-
    no.stdout(suppressWarnings(fitContinuous(phy, x, SE=se, model="lambda",
                                             control=list(niter=5))))

  ## Even though we don't agree on where the ML point is, let's check
  ## that the likelihood calculations do agree at the respective
  ## points.  First, at the diversitree ML point:
  p <- coef(fit.vcv)[c("lambda", "s2")]
  names(p)[2] <- "sigsq"
  expect_that(as.numeric(fit.geiger$lik(p)),
              equals(fit.vcv$lnLik))

  ## Second, at the geiger ML point.
  p <- coef(fit.geiger)[c("sigsq", "lambda")]
  expect_that(lik.vcv.se(p), equals(fit.geiger$opt$lnL))
  expect_that(lik.pru.R.se(p), equals(fit.geiger$opt$lnL))
  expect_that(lik.pru.C.se(p), equals(fit.geiger$opt$lnL))
})
