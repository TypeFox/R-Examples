context("Test rlassoEffect")
library(hdm)
library(testthat)

DGP_rlasso <- function(n, p, px){
  
  X <- matrix(rnorm(n*p), ncol=p)
  beta <- c(rep(2,px), rep(0,p-px))
  intercept <- 1
  y <- intercept + X %*% beta + rnorm(n)
  
  list(X = X, y = y, beta = beta)
}


set.seed(2)
ret <- DGP_rlasso(200, 100, 10)
X <- ret$X
y <- ret$y
beta <- ret$beta
frame <- as.data.frame(cbind(y, X))
colnames(frame) <- c("y", paste0("x", 1:100))
rm(ret)


test_that("rlassoEffect - Input check x, y and d",{
  expect_is(rlassoEffect(X, y, d = X[ ,1]), "rlassoEffects")
  expect_is(rlassoEffect(X, as.vector(y), d = X[ ,1]), "rlassoEffects")
  expect_is(rlassoEffect(X[, 1, drop = FALSE], y, d = X[ ,1]), "rlassoEffects")
  expect_is(rlassoEffect(X[, 1, drop = FALSE], as.vector(y), d = X[ ,1]), "rlassoEffects")
})

test_that("rlassoEffect - Input check I3",{
  expect_is(rlassoEffect(X, y, d = X[, 1], I3 = c(rep(TRUE,2),rep(FALSE,2), TRUE)), "rlassoEffects")
  expect_is(rlassoEffect(X, y, d = X[, 1], I3 = c(rep(TRUE,55),rep(FALSE,44), TRUE)), "rlassoEffects")
})

test_that("rlassoEffect - Input check post, intercept and normalize",{
  expect_is(rlassoEffect(X, y, d = X[ ,1], post = FALSE), "rlassoEffects")
  expect_is(rlassoEffect(X, y, d = X[ ,1], intercept = FALSE), "rlassoEffects")
  expect_is(rlassoEffect(X, y, d = X[ ,1], normalize = FALSE), "rlassoEffects")
  expect_is(rlassoEffect(X, y, d = X[ ,1], normalize = FALSE, intercept = FALSE), "rlassoEffects")
})

test_that("rlassoEffect - Input check penalty",{
  expect_is(rlassoEffect(X, y, d = X[ ,1], penalty = list(homoscedastic = FALSE, X.dependent.lambda =FALSE)), "rlassoEffects")
  expect_is(rlassoEffect(X, y, d = X[ ,1], penalty = list(homoscedastic = TRUE, X.dependent.lambda =FALSE)), "rlassoEffects")
  expect_is(rlassoEffect(X, y, d = X[ ,1], penalty = list(homoscedastic = FALSE, X.dependent.lambda =TRUE, numSim = 4000)), "rlassoEffects")
  expect_is(rlassoEffect(X, y, d = X[ ,1], penalty = list(homoscedastic = TRUE, X.dependent.lambda =TRUE, numSim = 4000)), "rlassoEffects")
  expect_is(rlassoEffect(X, y, d = X[ ,1], penalty = list(homoscedastic = "none", X.dependent.lambda =FALSE, lambda.start = 100)), "rlassoEffects")
  expect_is(rlassoEffect(X, y, d = X[ ,1], penalty = list(homoscedastic = "none", X.dependent.lambda =TRUE, lambda.start = 100)), "rlassoEffects")
  expect_is(rlassoEffect(X, y, d = X[ ,1], intercept = FALSE, penalty = list(homoscedastic = "none", X.dependent.lambda =FALSE, lambda.start = 100)), "rlassoEffects")
  expect_is(rlassoEffect(X, y, d = X[ ,1], penalty = list(homoscedastic = FALSE, X.dependent.lambda =FALSE, 
                                        lambda.start = NULL, c = 1.1, gamma = 0.1)), "rlassoEffects")
})

test_that("rlassoEffect - Input check control",{
  expect_is(rlassoEffect(X, y, d = X[ ,1], control = list(numIter = 15, tol = 10^-4, threshold = 10^-3)),"rlassoEffects")
  expect_is(rlassoEffect(X, y, d = X[ ,1], control = list(numIter = 25)), "rlassoEffects")         
})


