context("fastLasso - highdimensional regression")


## load packages
library("lars", quietly=TRUE)
library("robustHD", quietly=TRUE)

## set seed for reproducibility
set.seed(1234)

## generate data for tests
n <- 20                              # number of observations
p <- 50                              # number of predictors
beta <- rep.int(c(1, 0), c(5, p-5))  # true coefficients
x <- replicate(p, rnorm(n))          # predictor matrix
y <- c(x %*% beta + rnorm(n))        # response


## run tests

test_that("different values for lambda yield correct solution", {
  
  ## fit lasso with lars() as reference solution
  fitLars <- lars(x, y, type="lasso")
  
  ## extract values of lambda according to parametrization in robustHD
  lambda <- 2 * fitLars$lambda / n
  sMax <- length(lambda)
  
  ## choose different values of lambda and check solutions
  lambda <- c(lambda[1] * 1.5, 
              sort.int(union(lambda, (lambda[-sMax] + lambda[-1]) / 2), 
                       decreasing=TRUE), 
              lambda[sMax] * 0.5, 0.00001, 0)
  
  ## extract coefficients from solution computed via lars()
  coefLars <- sapply(n*lambda/2, function(l) {
    beta <- coef(fitLars, s=l, mode="lambda")
    alpha <- fitLars$mu - sum(beta * fitLars$meanx)
    c(alpha, beta)
  })
  
  ## fit models with fastLasso() and extract coefficients
  coefFastLasso <- sapply(lambda, function(l) {
    fitFastLasso <- robustHD:::fastLasso(x, y, lambda=l)
    coef(fitFastLasso)
  })
  
  ## test whether coefficients are equal
  expect_equal(coefLars, coefFastLasso)
})
