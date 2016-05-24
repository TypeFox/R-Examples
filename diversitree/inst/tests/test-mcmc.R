source("helper-diversitree.R")

context("MCMC")

make.mvn <- function(mean, vcv) {
  logdet <- as.numeric(determinant(vcv, TRUE)$modulus)
  tmp <- length(mean) * log(2 * pi) + logdet
  vcv.i <- solve(vcv)
  
  function(x) {
    dx <- x - mean
    -(tmp + rowSums((dx %*% vcv.i) * dx))/2
  }
}

vcv <- matrix(c(1, .25, .25, .75), 2, 2)
lik <- make.mvn(c(0, 0), vcv)

test_that("Basic mcmc works", {
  set.seed(1)
  n <- 100
  samples1 <- mcmc(lik, c(0, 0), n, 1, print.every=0)
  expect_that(nrow(samples1), equals(n))
  expect_that(names(samples1), is_identical_to(c("i", "X1", "X2", "p")))
})

test_that("Some options", {
  set.seed(1)
  expect_that(ign <- mcmc(lik, c(0, 0), 1, 1, print.every=1),
              prints_text("^[0-9]+: \\{[0-9,. -]+\\} -> [.0-9-]+$"))
})

test_that("Continuing a mcmc chain works", {
  set.seed(1)
  p0 <- c(0, 0)
  n1 <- 50
  n2 <- n1 + 50
  set.seed(1)  
  samples1 <- mcmc(lik, p0, n1, 1, print.every=0)
  samples2 <- mcmc(lik, NULL, n2, 1, print.every=0,
                   previous=samples1)

  ## We sample *up* to n2, not for n2 more points.
  expect_that(nrow(samples2), equals(n2))

  ## Restarting is deterministic if the RNG is in the right place.
  set.seed(1)
  samples3 <- mcmc(lik, p0, n2, 1, print.every=0)
  expect_that(samples3, is_identical_to(samples2))

  ## Can't provide both a starting point and previous samples.
  expect_that(mcmc(lik, c(0,0), n2, 1, print.every=0,
                   previous=samples1), throws_error())
})

test_that("Likelihood function is saved with fit", {
  samples <- mcmc(lik, c(0, 0), 10, 1, print.every=0)
  expect_that(samples, has_attribute("func"))
  expect_that(attr(samples, "func"), is_a("function"))
  expect_that(attr(samples, "func"), equals(lik))

  samples.no.func <- mcmc(lik, c(0, 0), 10, 1, print.every=0,
                          keep.func=FALSE)
  expect_that(attr(samples.no.func, "func"), is_null())

  expect_that(attr(drop.likelihood(samples),         "func"), is_null())
  expect_that(attr(drop.likelihood(samples.no.func), "func"), is_null())
})

test_that("dtlik mcmcsamples have correct class", {
  samples <- mcmc(lik, c(0, 0), 10, w=1, print.every=0)
  expect_that(class(samples),
              equals(c("mcmcsamples", "data.frame")))
  
  pars <- c(0.1, 0.03)
  set.seed(2)
  phy <- tree.bd(pars, max.taxa=60)
  lik <- make.bd(phy)
  samples <- mcmc(lik, pars, 10, w=1, print.every=0)
  expect_that(class(samples),
              equals(c("mcmcsamples.bd", "mcmcsamples", "data.frame")))
})

test_that("Argument modification is saved at function save", {
  pars <- c(0.1, 0.03)
  set.seed(2)
  phy <- tree.bd(pars, max.taxa=60)
  lik <- make.bd(phy)
  ## Otherwise the stuff below has no effect.
  expect_that(formals(lik)$condition.surv, is_true())

  samples <- mcmc(lik, pars, 10, w=1, print.every=0,
                  condition.surv=FALSE)
  expect_that(samples, has_attribute("func"))
  expect_that(formals(attr(samples, "func"))$condition.surv,
              is_false())

  expect_that(samples, has_attribute("func"))
  expect_that(formals(attr(samples, "func"))$condition.surv,
              is_false())
  ## Will be simplified by the new "devtools::not()".
  expect_that(identical(attr(samples, "func"), lik), is_false())
})

## Not tested:
## 
##   * Correct behaviour when sampling from univariate models
##     (i.e. check that drop=FALSE works)
##   * Check that expansion of constrained models works via
##     full=TRUE.
test_that("coef.mcmcsamples works", {
  n <- 300
  set.seed(1)
  samples <- mcmc(lik, c(0, 0), n, w=5, print.every=0)
  p <- coef(samples)

  expect_that(samples, is_a("data.frame"))
  expect_that(p, is_a("matrix"))
  expect_that(ncol(p), equals(ncol(samples) - 2))
  expect_that(nrow(p), equals(n))

  nb <- 100
  p <- coef(samples, burnin=nb)
  expect_that(ncol(p), equals(ncol(samples) - 2))
  expect_that(nrow(p), equals(n - nb))

  pb <- 1/pi
  p <- coef(samples, burnin=pb)
  expect_that(nrow(p), equals(n - floor(n * pb)))
  
  thin <- 7
  p <- coef(samples, thin=thin)
  expect_that(ncol(p), equals(ncol(samples) - 2))
  expect_that(nrow(p), equals(ceiling(n / thin)))

  p <- coef(samples, burnin=nb, thin=thin)
  expect_that(ncol(p), equals(ncol(samples) - 2))
  expect_that(nrow(p), equals(ceiling((n - nb) / thin)))

  ns <- 50
  p <- coef(samples, sample=ns)
  expect_that(ncol(p), equals(ncol(samples) - 2))
  expect_that(nrow(p), equals(ns))

  ns2 <- floor((n - nb) / thin)
  expect_that(p <- coef(samples, burnin=nb, thin=thin, sample=ns),
              gives_warning())
  
  p <- coef(samples, burnin=nb, thin=thin, sample=ns2)
  expect_that(nrow(p), equals(ns2))  
})

## Not tested yet:

## Then, with saving:
## samples1 <- mcmc(lik, c(0, 0), 10000, 1, print.every=1000,
##                  save.every=10, save.file="test.rds")

## samples1 <- mcmc(lik, c(0, 0), 10000, 1, print.every=1000,
##                  save.every.dt=seconds(4), save.file="test.rds")

## file.remove("test.rds")
