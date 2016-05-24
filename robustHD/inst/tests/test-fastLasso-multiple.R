context("fastLasso - multiple regression")


## load packages
library("lars", quietly=TRUE)
library("robustHD", quietly=TRUE)

## set seed for reproducibility
set.seed(1234)

## generate data for tests
n <- 20                        # number of observations
beta <- c(1, 1, 0, 0, 0)       # true coefficients
p <- length(beta)              # number of predictors
x <- replicate(p, rnorm(n))    # predictor matrix
y <- c(x %*% beta + rnorm(n))  # response


## run tests

test_that("special case for no penalty yields LS solution", {
  
  ## compute LS solution and extract coefficients
  fitLS <- lm(y~x)
  coefLS <- unname(coef(fitLS))
  
  ## fit models with fastLasso() and extract coefficients
  fitFastLasso <- robustHD:::fastLasso(x, y, lambda=0)
  coefFastLasso <- coef(fitFastLasso)
  
  ## test whether coefficients are equal
  expect_equal(coefLS, coefFastLasso)
})


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
              lambda[sMax] * 0.5, 0.00001)
  
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
