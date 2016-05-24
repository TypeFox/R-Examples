context("createTransform")

# check what the transform does to the zero and unit vectors
invariants <- function(basis) {
  n <- nrow(basis$basis)
  m <- ncol(basis$basis)
  zero <- rep(0, m)
  unit <- rep(sqrt(1/m), m)

  transform <- createTransform(basis, inverse=FALSE, keepHomogeneous=FALSE)
  expect_equal(
    as.vector(transform %*% c(zero, 1)),
    basis$translate)
  expect_equal(
    as.vector(transform %*% c(unit, 1)),
    as.vector(basis$basis %*% unit + basis$translate))

  transform <- createTransform(basis, inverse=FALSE, keepHomogeneous=TRUE)
  expect_equal(
    as.vector(transform %*% c(zero, 1)),
    c(basis$translate, 1))
  expect_equal(
    as.vector(transform %*% c(unit, 1)),
    c(basis$basis %*% unit + basis$translate, 1))

  inverse <- createTransform(basis, inverse=TRUE, keepHomogeneous=FALSE)
  expect_equal(
    as.vector(inverse %*% c(basis$translate, 1)),
    zero)
  expect_equal(
    as.vector(inverse %*% c(basis$basis %*% unit + basis$translate, 1)),
    unit)

  inverse <- createTransform(basis, inverse=TRUE, keepHomogeneous=TRUE)
  expect_equal(
    as.vector(inverse %*% c(basis$translate, 1)),
    c(zero, 1))
  expect_equal(
    as.vector(inverse %*% c(basis$basis %*% unit + basis$translate, 1)),
    c(unit, 1))
}

n <- 3

test_that("identity transform", {
  basis <- list(basis=diag(n), translate=rep(0, n))
  invariants(basis)
  transform <- createTransform(basis, inverse=FALSE, keepHomogeneous=TRUE)
  expect_equal(transform, diag(n+1))
  transform <- createTransform(basis, inverse=TRUE, keepHomogeneous=TRUE)
  expect_equal(transform, diag(n+1))

  transform <- createTransform(basis, inverse=FALSE, keepHomogeneous=FALSE)
  expect_equal(transform, cbind(diag(n), rep(0, n)))
  transform <- createTransform(basis, inverse=TRUE, keepHomogeneous=FALSE)
  expect_equal(transform, cbind(diag(n), rep(0, n)))
})

test_that("translation only", {
  basis <- list(basis=diag(n), translate=rep(1/2, n))
  invariants(basis)
  transform <- createTransform(basis, inverse=FALSE, keepHomogeneous=TRUE)
  expected <- diag(n+1)
  expected[1:3,4] <- 1/2
  expect_equal(transform, expected)
  transform <- createTransform(basis, inverse=FALSE, keepHomogeneous=FALSE)
  expect_equal(transform, expected[1:3,])
  transform <- createTransform(basis, inverse=TRUE, keepHomogeneous=TRUE)
  expected[1:3,4] <- -1/2
  transform <- createTransform(basis, inverse=TRUE, keepHomogeneous=FALSE)
  expect_equal(transform, expected[1:3,])
})

test_that("reduce dimensionality by 1", {
  basis <- list(basis=cbind(c(1,0,0),sqrt(1/2)*c(0,-1,1)), translate=rep(0,n))
  invariants(basis)
  transform <- createTransform(basis, inverse=FALSE, keepHomogeneous=TRUE)
  expected <- matrix(0, nrow=4, ncol=3)
  expected[1:3,1:2] <- basis$basis
  expected[4,3] <- 1
  expect_equal(transform, expected)
  transform <- createTransform(basis, inverse=FALSE, keepHomogeneous=FALSE)
  expect_equal(transform, expected[1:3,])

  transform <- createTransform(basis, inverse=TRUE, keepHomogeneous=TRUE)
  expected <- matrix(0, nrow=3, ncol=4)
  expected[1:2,1:3] <- t(basis$basis)
  expected[3,4] <- 1
  expect_equal(transform, expected)
  transform <- createTransform(basis, inverse=TRUE, keepHomogeneous=FALSE)
  expect_equal(transform, expected[1:2,])
})

test_that("reduce dimensionality by 2", {
  basis <- list(basis=cbind(sqrt(1/3)*c(1,-1,1)), translate=rep(0,n))
  invariants(basis)
  transform <- createTransform(basis, inverse=FALSE, keepHomogeneous=TRUE)
  expected <- matrix(0, nrow=4, ncol=2)
  expected[1:3,1] <- basis$basis
  expected[4,2] <- 1
  expect_equal(transform, expected)
  transform <- createTransform(basis, inverse=FALSE, keepHomogeneous=FALSE)
  expect_equal(transform, expected[1:3,])

  # inverse failed in <= 0.4-1
  transform <- createTransform(basis, inverse=TRUE, keepHomogeneous=TRUE)
  expected <- matrix(0, nrow=2, ncol=4)
  expected[1,1:3] <- t(basis$basis)
  expected[2,4] <- 1
  expect_equal(transform, expected)
  transform <- createTransform(basis, inverse=TRUE, keepHomogeneous=FALSE)
  expect_equal(transform, expected[1,,drop=FALSE])
})

test_that("reduce dimensionality and translate", {
  basis <- list(basis=cbind(sqrt(1/3)*c(1,-1,1)), translate=c(1/2,1/3,1/6))
  invariants(basis)
  transform <- createTransform(basis, inverse=FALSE, keepHomogeneous=TRUE)
  expected <- matrix(0, nrow=4, ncol=2)
  expected[1:3,1] <- basis$basis
  expected[1:3,2] <- basis$translate
  expected[4,2] <- 1
  expect_equal(transform, expected)
  transform <- createTransform(basis, inverse=FALSE, keepHomogeneous=FALSE)
  expect_equal(transform, expected[1:3,])

  transform <- createTransform(basis, inverse=TRUE, keepHomogeneous=TRUE)
  expected <- matrix(0, nrow=2, ncol=4)
  expected[1,1:3] <- t(basis$basis)
  expected[1,4] <- -mean(basis$translate) * sqrt(1/3)
  expected[2,4] <- 1
  expect_equal(transform, expected)
  transform <- createTransform(basis, inverse=TRUE, keepHomogeneous=FALSE)
  expect_equal(transform, expected[1,,drop=FALSE])
})
