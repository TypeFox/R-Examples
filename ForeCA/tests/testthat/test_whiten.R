
kNumVariables <- 3
kNumObs <- 100
kSeries <- matrix(rnorm(kNumObs * kNumVariables), ncol = kNumVariables)
kSeriesCentered <- sweep(kSeries, 2, colMeans(kSeries), "-")
kPosDef <- cov(kSeriesCentered)

whitening.result <- whiten(kSeries)
whitening.result.centered <- whiten(kSeriesCentered)

test_that("whiten has the right names", {
  # check names of the whitening.result list
  expect_equal(names(whitening.result), 
               c("center", "scale", "values", "whitening", "dewhitening", "U"))
  expect_equal(names(whiten(whitening.result$U)), 
               c("center", "scale", "values", "whitening", "dewhitening", "U"))
  expect_true(attr(whitening.result$U, "whitened"),
              info = "whitened attributed is TRUE")
})

test_that("whiten does indeed give identity covariance matrix", {
  # covariance is the identity matrix and mean is 0
  expect_equal(cov(whitening.result$U), diag(1, kNumVariables))
  expect_equal(colMeans(whitening.result$U), rep(0, ncol(kSeries)))
})

test_that("whiten gives same normalized data from centered vs uncentered input", {
  expect_equal(whitening.result.centered$U, whitening.result$U)
})

test_that("whitening transformation and inverse are indeed inverse transformations", {
  expect_equal(whitening.result$whitening %*% whitening.result$dewhitening,  
               diag(1, kNumVariables))
})

test_that("dewhitening transformation rows give the correct column from original data", {
  expect_equal(cor(c(whitening.result$U %*% whitening.result$dewhitening[1, ]),
                   kSeries[, 1]), 1)
})

test_that("dewhitening transformation rows give the correct column from original data", {
  expect_equal(cor(c(whitening.result$U %*% whitening.result$dewhitening),
                   c(scale(kSeries, scale = FALSE))), 1)
})

test_that("whitening transformation column gives the correct column from U", {
  expect_equal(cor(c(kSeries %*% whitening.result$whitening[, 1]),
                   whitening.result$U[, 1]), 1)
})

context("whiten for univariate series")
one.var.whiten <- whiten(kSeries[, 1])
test_that("mean is 0", {
  expect_equal(0, mean(one.var.whiten$U))
})

test_that("var is 1", {
  expect_equal(1, var(one.var.whiten$U))
})

test_that("whitening is inverse of dewhitening", {
  expect_equal(one.var.whiten$whitening, 1 / one.var.whiten$dewhitening)
})

context("Testing check_whitened function")

test_that("check_whitened returns error for non-zero mean input", {
  # must have 0 mean
  expect_error(check_whitened(kSeries, FALSE))
})

test_that("check_whitened returns error for not unit variance input", {
  # must have unit variance
  expect_error(check_whitened(scale(kSeries, scale = FALSE), FALSE))
})

test_that("check_whitened returns error for correlated input", {
  # must have 0 mean
  expect_error(check_whitened(scale(kSeries, scale = TRUE), FALSE))
})

test_that("check_whitened returns does nothing for whitened input", {
  expect_silent(check_whitened(whiten(kSeries)$U))
  expect_silent(check_whitened(whiten(kSeries)$U, FALSE))
})

test_that("check_whitened also works for univariate input", {
  univ.u <- whiten(kSeries)$U[, 1]
  expect_silent(check_whitened(univ.u, FALSE))
  attr(univ.u, "whitened") <- TRUE
  expect_silent(check_whitened(univ.u, TRUE))
})


context("Testing sqrt_matrix function")

test_that("sqrt_matrix takes only square matrix input", {
  # must be a matrix
  expect_error(sqrt_matrix(c(1, 2, 3)))
  # must be square
  expect_error(sqrt_matrix(kSeriesCentered))
})

test_that("sqrt_matrix(A) computes in fact the square root B*B = A", {
  sqrt.mat <- sqrt_matrix(kPosDef)
  
  # it is actually the sqaare root
  expect_equal(sqrt.mat %*% sqrt.mat, kPosDef)
  # diagonal is its own square root
  expect_equal(diag(1, 10), sqrt_matrix(diag(1, 10)))
})

test_that("sqrt_matrix handles negative eigenvalues / complex values correctly", {
  # complex values
  expect_equal(diag(0 + (0+1i), 10), sqrt_matrix(diag(-1, 10)))
})

context("sqrt_matrix auxiliary metrics (transformation/eigenvalues")
sqrt.result <- sqrt_matrix(kPosDef, return.sqrt.only = FALSE, 
                           symmetric = TRUE)

test_that("sqrt_matrix returns a list with correct names", {
  # return more than just the inverse
  expect_true(inherits(sqrt.result, "list"))
  expect_equal(names(sqrt.result), c("values", "vectors", "sqrt", "sqrt.inverse"))
})

test_that("sqrt_matrix gives correct inverse for whitening", {
  # sqrt.inverse is actually the inverse to sqrt
  expect_equal(sqrt.result$sqrt %*%  sqrt.result$sqrt.inverse, 
               diag(1, ncol(kPosDef)))
  
  expect_equal(kPosDef %*% sqrt.result$sqrt.inverse %*%  sqrt.result$sqrt.inverse,        
               diag(1, ncol(kPosDef)))
})


