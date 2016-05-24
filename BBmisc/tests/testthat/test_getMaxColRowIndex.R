context("getMaxIndexOfRows")

test_that("getMaxIndexOfRows", {
  a = matrix(1:6, nrow=2)
  expect_equal(getMaxIndexOfRows(a), c(3L, 3L))
  a = matrix(6:1, nrow=2)
  expect_equal(getMaxIndexOfRows(a), c(1L, 1L))
  a = rbind(c(1, 999), c(-1, -5))
  expect_equal(getMaxIndexOfRows(a), c(2L, 1L))
  a = matrix(rnorm(50*10), nrow=50)
  expect_equal(getMaxIndexOfRows(a), apply(a, 1, which.max))
})

test_that("getMaxIndexOfCols", {
  a = matrix(1:6, nrow=2)
  expect_equal(getMaxIndexOfCols(a), c(2L, 2L, 2L))
  a = matrix(6:1, nrow=2)
  expect_equal(getMaxIndexOfCols(a), c(1L, 1L, 1L))
  a = rbind(c(1, 999), c(-1, -5))
  expect_equal(getMaxIndexOfCols(a), c(1L, 1L))
  a = matrix(rnorm(50*10), nrow=50)
  expect_equal(getMaxIndexOfCols(a), apply(a, 2, which.max))
})

test_that("normal", {
  expect_equal(getMaxIndexOfRows(diag(10)), 1:10)
  n = 100
  perm = sample(n)
  D = diag(n)
  expect_equal(getMaxIndexOfRows(D[perm, ]), (1:n)[perm])
})

test_that("NA values", {
  n = 300
  m = matrix(runif(n), ncol=3)
  mm = m
  mm[, 2] = NA
  expect_equal(getMaxIndexOfRows(mm), rep(NA_integer_, n/3))
  
  a = matrix(c(1, NA, 2, 3, NA, NA), nrow=3, byrow=TRUE)
  expect_equal(getMaxIndexOfRows(a, na.rm=FALSE), c(NA, 2L, NA))  
  expect_equal(getMaxIndexOfRows(a, na.rm=TRUE), c(1L, 2L, -1))  
})

test_that("infinite values", {
  n = 300
  m = matrix(runif(n), ncol=3)
  m[, 2] = Inf
  expect_equal(getMaxIndexOfRows(m), rep(2L, 100L))
})

test_that("max.col oddity", {
  expect_equal(getMaxIndexOfRows(cbind(1:10, 2:11, -Inf)), rep(2, 10))
  expect_equal(getMaxIndexOfRows(cbind(-1e9 * 1:10, 1:10, 2:11)), rep(3, 10))
})

test_that("ties", {
  a = matrix(c(1, 1, 2, 2), nrow=2, byrow=TRUE)
  expect_equal(getMaxIndexOfRows(a, ties.method="first"), c(1L, 1L))
  expect_equal(getMaxIndexOfRows(a, ties.method="last"), c(2L, 2L))
  a = matrix(c(2, 1, 2, 2, 2, 1), nrow=2, byrow=TRUE)
  expect_equal(getMaxIndexOfRows(a, ties.method="first"), c(1L, 1L))
  expect_equal(getMaxIndexOfRows(a, ties.method="last"), c(3L, 2L))
  a = matrix(c(1, 1, 2, 2), nrow=2, byrow=TRUE)
  expect_equal(getMaxIndexOfCols(a, ties.method="first"), c(2L, 2L))
  expect_equal(getMaxIndexOfCols(a, ties.method="last"), c(2L, 2L))
  a = matrix(c(2, 1, 2, 2, 2, 1), nrow=2, byrow=TRUE)
  expect_equal(getMaxIndexOfCols(a, ties.method="first"), c(1L, 2L, 1L))
  expect_equal(getMaxIndexOfCols(a, ties.method="last"), c(2L, 2L, 1L))
})


