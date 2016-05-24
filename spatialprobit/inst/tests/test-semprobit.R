context("semprobit")

library(spatialprobit)
set.seed(2)
n <- 200
I_n <- sparseMatrix(i=1:n, j=1:n, x=1)
rho <- 0.75
W <- kNearestNeighbors(x=rnorm(n), y=rnorm(n))
rho <- 0.5
S <- (I_n - rho * W)
QR <- qr(S)
x <- rnorm(n)

test_that("semprobit() only accepts valid inputs for y", {
  # y must only have 0 or 1 values
  y <- rep(0, n)
  y[1] <- 2
  expect_that(semprobit(y ~ x, W=W), throws_error())
  
  # y must have same length as x
  y <- rep(0, n-1)
  expect_that(semprobit(y ~ x, W=W), throws_error())
})

test_that("semprobit() only accepts valid inputs for W", {
  
  # only dummy values for y as we expect errors to be thrown
  y <- rep(0, n)
    
  # spatial weights matrix W must be a matrix
  W <- "A"
  expect_that(semprobit(y ~ x, W=W), throws_error())
  
  # spatial weights matrix must be a sparse matrix
  W <- diag(n)
  expect_that(semprobit(y ~ x, W=W), throws_error())
  
  # spatial weights matrix W must not contain non-zeros in the diagonal
  W <- Matrix(diag(n), sparse=TRUE)
  expect_that(semprobit(y ~ x, W=W), throws_error())
  
  W <- matrix(0, n, n)
  W[1,1] <- 1
  expect_that(semprobit(y ~ x, W=W), throws_error())
})