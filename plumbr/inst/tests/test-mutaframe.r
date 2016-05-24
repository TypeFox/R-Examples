library(testthat)

context("Basic construction and mutability")

test_that("can construct empty mutaframe", {
  empty <- mutaframe()
  
  expect_that(nrow(empty), equals(0))
  expect_that(ncol(empty), equals(0))
})

test_that("can construct single column mutaframes", {
  a <- mutaframe(1:10)
  expect_that(ncol(a), equals(1))
  
  b <- mutaframe(a = 1:10)
  expect_that(ncol(a), equals(1))
})

test_that("shorter columns recycled", {
  a <- mutaframe(a = 1:10, b = 1:5)
  expect_that(nrow(a), equals(10))

  b <- mutaframe(a = 1:10, b = 1:5, c = 1:2)
  expect_that(nrow(b), equals(10))
})

test_that("failure to recycle raises error", {
  expect_that(mutaframe(a = 1:10, b = 1:9), throws_error("different row counts"))
})

test_that("default names behave like data.frame", {
  a <- mutaframe(1, 2)
  expect_that(names(a), equals(c("X1", "X2")))
  
  b <- mutaframe(1, X1 = 2)
  expect_that(names(b), equals(c("X1", "X1.1")))

  c <- mutaframe(`a b` = 1)
  expect_that(names(c), equals("a.b"))
  
})

test_that("can construct from a dataframe", {
  m <- mutaframe(mtcars)

  expect_that(names(m), equals(names(mtcars)))
  expect_that(dim(m), equals(dim(mtcars)))
  
  expect_that(m[, 1], equals(mtcars[, 1]))
  expect_that(m[, 5], equals(mtcars[, 5]))
})

test_that("coercing data frame works", {
  m <- as.mutaframe(mtcars)

  expect_that(names(m), equals(names(mtcars)))
  expect_that(dim(m), equals(dim(mtcars)))
  
  expect_that(m[, 1], equals(mtcars[, 1]))
  expect_that(m[, 5], equals(mtcars[, 5]))
})