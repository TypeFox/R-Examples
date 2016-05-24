source("helper-diversitree.R")

context("MLE")

## Will be useful to have a function for testing tolerance to within
## 1e-7, as that works out to be how accurate most things actually
## are.

## Start with a simple 2 parameter model from the BD model
pars <- c(0.1, 0.03)
set.seed(2)
phy <- tree.bd(pars, max.taxa=60)
lik <- make.bd(phy)

## First, exercise some basic options:
ans.nlm.1 <- find.mle(lik, pars, method="nlm")
expect_that(ans.nlm.2 <- find.mle(lik, pars, method="nlm",
                                  verbose=40),
            prints_text("."))
ans.nlm.3 <- find.mle(lik, pars, method="nlm", fail.value=-500)
ans.nlm.4 <- find.mle(lik, pars, method="nlm",
                      control=list(fail.penalty=500))

expect_that(ans.nlm.1, is_identical_to(ans.nlm.2))
expect_that(ans.nlm.1[1:2], equals(ans.nlm.3[1:2], tolerance=1e-2))
expect_that(!identical(ans.nlm.1, ans.nlm.4), is_true())
expect_that(ans.nlm.1[1:2], equals(ans.nlm.4[1:2], tolerance=0.07))

## Next, excercise the different methods:
ans.optim <- suppressWarnings(find.mle(lik, pars, method="optim"))
ans.subplex <- find.mle(lik, pars, method="subplex")
ans.nlminb <- find.mle(lik, pars, method="nlminb")
ans.nlm <- find.mle(lik, pars, method="nlm")
ans.minqa <- find.mle(lik, pars, method="minqa")

## Different optim methods:
ans.optim.nm <- find.mle(lik, pars, method="optim",
                         control=list(optim.method="Nelder-Mead"))
ans.optim.bfgs <- find.mle(lik, pars, method="optim",
                           control=list(optim.method="BFGS"))
ans.optim.cg <- find.mle(lik, pars, method="optim",
                         control=list(optim.method="CG"))
ans.optim.lbfgsb <-
  suppressWarnings(find.mle(lik, pars, method="optim",
                            control=list(optim.method="L-BFGS-B")))
## Just takes too long, and not interesting.
## ans.optim.sann <- find.mle(lik, pars, method="optim",
##                            control=list(optim.method="SANN"))

## Different minqa methods:
ans.minqa.n <- find.mle(lik, pars, method="minqa",
                        control=list(minqa.method="newuoa"))
ans.minqa.b <- find.mle(lik, pars, method="minqa",
                        control=list(minqa.method="bobyqa"))
ans.minqa.u <- find.mle(lik, pars, method="minqa",
                        control=list(minqa.method="uobyqa"))

equals.tol <- function(expected, tol, ...)
  equals(expected, tolerance=tol, ...)

tol <- 0.11
ans <- c(lambda=0.0574260691275107, mu=0.0)

expect_that(coef(ans.optim), equals.tol(ans, tol))
expect_that(coef(ans.subplex), equals.tol(ans, tol))
expect_that(coef(ans.nlminb), equals.tol(ans, tol))
expect_that(coef(ans.nlm), equals.tol(ans, tol))
expect_that(coef(ans.minqa), equals.tol(ans, tol))

expect_that(coef(ans.optim.nm), equals.tol(ans, tol))
expect_that(coef(ans.optim.bfgs), equals.tol(ans, tol))
## Looks like CG fails on this problem.
## expect_that(coef(ans.optim.cg), equals.tol(ans, tol))
expect_that(coef(ans.optim.lbfgsb), equals.tol(ans, tol))
## expect_that(coef(ans.optim.sann), equals.tol(ans, tol))

expect_that(coef(ans.minqa.n), equals.tol(ans, tol))
expect_that(coef(ans.minqa.b), equals.tol(ans, tol))
expect_that(coef(ans.minqa.u), equals.tol(ans, tol))

## "Odd" methods, involving integers.

## Here is a "likelihood" function based on the Rosenbrock banana
## function.  We will evaluate this function only at integer values of
## x.
lik <- function(x, x2, as.integer=TRUE) {
  x1 <- if ( as.integer ) diversitree:::check.integer(x) else x
  -(100*(x2-x1*x1)^2+(1-x1)^2)
}

fit <- find.mle(lik, -11, x2=-33, method="int1d")
expect_that(fit$par, is_identical_to(0))

## And a few options here, too:
fit.1 <- find.mle(lik, -11, x2=-33, method="int1d", upper=11)
fit.2 <- find.mle(lik, -11, x2=-33, method="int1d", lower=-11)
fit.3 <- find.mle(lik, -11, x2=-33, method="int1d",
                  control=list(interval=c(-11, 1)))
expect_that(fit.1$par, is_identical_to(0))
expect_that(fit.2$par, is_identical_to(0))
expect_that(fit.3$par, is_identical_to(0))

## Now, the mixed method: the 4d rosenbrock function with some integer
## axes.
rosen.multi <- function(x) {
  xx <- matrix(x, 2)
  x1 <- xx[1,]
  x2 <- xx[2,]
  -sum(100*(x2-x1*x1)^2+(1-x1)^2)
}

## Let's have the third argument be an integer.
rosen.mixed <- function(x) {
  diversitree:::check.integer(x[3])
  rosen.multi(x)
}

## This works well.
p <- rep(-11, 4)
fit.c <- find.mle(rosen.multi, p, method="subplex")

## However, this will fail to fit, as the integer axis is highly
## correlated with one real axis (in fact, this confuses subplex, and
## the other rosenbrock function also does not complete).
fit.m <- find.mle(rosen.mixed, p, method="mixed",
                  control=list(is.integer=3))
fit.m2 <- find.mle(rosen.mixed, fit.m$par, method="mixed",
                   control=list(is.integer=3))
expect_that(fit.m2$par[3], is_identical_to(-11))

## Options around saving the objective function.
test_that("Likelihood function is saved with fit", { 
  lik <- make.bd(phy)
  pars <- c(0.1, 0.03)  
  fit <- find.mle(lik, pars)
  expect_that(fit, has_attribute("func"))
  expect_that(attr(fit, "func"), is_a("function"))
  expect_that(attr(fit, "func"), is_identical_to(lik))

  fit.no.func <- find.mle(lik, pars, keep.func=FALSE)
  expect_that(attr(fit.no.func, "func"), is_null())
  expect_that(get.likelihood(fit.no.func), throws_error())

  expect_that(attr(drop.likelihood(fit),         "func"), is_null())
  expect_that(attr(drop.likelihood(fit.no.func), "func"), is_null())

  expect_that(get.likelihood(drop.likelihood(fit)),
              throws_error())
})

test_that("Argument modification is saved at function save", { 
  lik <- make.bd(phy)
  ## Otherwise the stuff below has no effect.
  expect_that(formals(lik)$condition.surv, is_true())
  
  pars <- c(0.1, 0.03)  
  fit <- find.mle(lik, pars, condition.surv=FALSE)
  expect_that(fit, has_attribute("func"))
  expect_that(formals(attr(fit, "func"))$condition.surv,
              is_false())

  ## Will be simplified by the new "devtools::not()".
  expect_that(identical(attr(fit, "func"), lik), is_false())
})
