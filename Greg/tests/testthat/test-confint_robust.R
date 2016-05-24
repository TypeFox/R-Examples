library(testthat)
library("sandwich")

context("Check confint_robust")
test_that("Check regular lm", {
  set.seed(100)  
  n <- 50
  x <- runif(n)
  y <- x + rnorm(n)
  
  fit <- lm(y~x)
  expect_equivalent(dim(confint_robust(fit, HC_type = "HC4m")), c(2,2))
})


test_that("Check robcov_alt", {
  set.seed(120)  
  n <- 50
  x <- runif(n)
  y <- x + rnorm(n)
  
  fit <- lm(y~x)
  expect_error(robcov_alt(fit, cluster = 2))
  
  expect_equal(dim(robcov_alt(fit)$var),
               dim(vcov(fit, intercepts = "all")))
})

test_that("Check hatvalues for ols", {
  set.seed(120)  
  n <- 50
  x <- runif(n)
  y <- x + rnorm(n)
  
  fit <- ols(y~x)
  expect_equivalent(hatvalues(fit),
                    ols.influence(fit)$hat)
})

test_that("Check bread for ols", {
  set.seed(120)  
  n <- 50
  x <- runif(n)
  y <- x + rnorm(n)
  
  fit <- ols(y~x, x = TRUE)
  expect_equal(dim(bread(fit)), c(2,2))
})

test_that("Check estfun for ols", {
  set.seed(120)  
  n <- 50
  x <- runif(n)
  y <- x + rnorm(n)
  
  fit <- ols(y~x, x = TRUE)
  expect_equal(nrow(estfun(fit)), 50)
})

test_that("Check model.matrix for ols", {
  set.seed(120)  
  n <- 50
  x <- runif(n)
  y <- x + rnorm(n)
  
  fit <- ols(y~x, x = TRUE)
  expect_equivalent(model.matrix(fit), 
                    cbind(rep(1, n),
                          x))

  fit <- ols(y~x)
  expect_warning(model.matrix(fit))
})
