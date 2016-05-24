context("Test miscellanous")

test_that("Beta integration range", {
  f <- function(x) 1/x^.2
  c <- 0.5; d<- 100
  range <- brr:::beta_integration_range(c,d,f) 
  I1 <- integrate(function(x) f(x)*dbeta(x,c,d), lower=0, upper=1)
  I2 <- integrate(function(x) f(x)*dbeta(x,c,d), lower=0, upper=range[2])
  expect_equal(I1$value, I2$value, tolerance=1e-5)
  #
  g <- function(x) f(1-x)
  c <- 100; d<- 0.5
  range <- brr:::beta_integration_range(c,d,g) 
  I3 = integrate(function(x) g(x)*dbeta(x,c,d), lower=range[1], upper=1)
  expect_equal(I2$value, I3$value, tolerance=1e-9)  
  # 
  c <- 1
  expect_identical(brr:::beta_integration_range(c,d,f), c(0,1))
  #
  c <- 100; d <- 10
  f <- function(x) dbeta(x, 0.5, 0.5)
  range <- brr:::beta_integration_range(c,d,f) 
  I1 <- integrate(function(x) f(x)*dbeta(x,c,d), lower=0, upper=1)
  I2 <- integrate(function(x) f(x)*dbeta(x,c,d), lower=0, upper=range[2])
  expect_equal(I1$value, I2$value, tolerance=1e-5)
  #
  f <- function(x) dbeta(x, 15, 15)
  range <- brr:::beta_integration_range(c,d,f) 
  I1 <- integrate(function(x) f(x)*dbeta(x,c,d), lower=0, upper=1)
  I2 <- integrate(function(x) f(x)*dbeta(x,c,d), lower=0, upper=range[2])
  expect_equal(I1$value, I2$value, tolerance=1e-5)  
})

test_that("dd_moment", {
  set.seed(666)
  tol <- 1e-10
  # Poisson 
  lambda <- rnorm(3, 5)
  m <- sapply(lambda, function(lambda) dd_moment(dpois, lambda=lambda, accuracy=tol))
  expect_equal(m, lambda, tolerance=tol)
  m2 <- sapply(lambda, function(lambda) dd_moment(dpois, k=2, lambda=lambda, accuracy=tol))
  expect_equal(m2-m^2, lambda, tolerance=tol)
  # beta_nbinom
  accuracy <- 1e-12
  tol <- 1e-9
  a <- rnorm(3, 5)
  c <- 6; d <- 3
  m <- sapply(a, function(a) dd_moment(dbeta_nbinom, a=a, c=c, d=d, accuracy=accuracy))
  m2 <- sapply(a, function(a) dd_moment(dbeta_nbinom, k=2, a=a, c=c, d=d, accuracy=accuracy))
  m1 <- sapply(a, function(a) summary_beta_nbinom(a=a, c=c, d=d)$mean) 
  SD <- sapply(a, function(a) summary_beta_nbinom(a=a, c=c, d=d)$sd) 
  expect_equal(m, m1, tolerance=tol)
  expect_equal(m2-m^2, SD^2, tolerance=tol)
})