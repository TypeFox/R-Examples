library(testthat)

context("Accessors")

test_that("getting and setting a single value works", {
  a <- mutaframe(a = 1:10)
  expect_that(a$a[1], equals(1))
  expect_that(a[[1]][1], equals(1))
  expect_that(a[["a"]][1], equals(1))
  expect_that(a[1][1 ,1], equals(1))
  expect_that(a["a"][1 ,1], equals(1))
  expect_that(a[1,1], equals(1))
  
  a[1, 1] <- 2
  expect_that(a[1,1], equals(2))
})


test_that("you can add a new column", {
  a <- mutaframe(a = 1:10)
  a$b <- 1:10 * 2
  expect_that(ncol(a), equals(2))
  expect_that(names(a), equals(c("a", "b")))
  expect_that(a$b, equals(1:10 * 2))
  
  a <- mutaframe(a = 1:10)
  a[, "b"] <- 1:10 * 2
  expect_that(ncol(a), equals(2))
  expect_that(names(a), equals(c("a", "b")))
  expect_that(a$b, equals(1:10 * 2))
})

test_that("you can delete a column", {
  a <- mutaframe(a = 1:10, b = 1:10)
  expect_that(ncol(a), equals(2))
  a$b <- NULL
  expect_that(ncol(a), equals(1))
  expect_that(names(a), equals("a"))
})
