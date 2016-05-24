context("is.empty")

test_that("non-empty vector", {
  a <- c(1,2,3)
  expect_that(is.empty(a), equals(FALSE))
})
test_that("empty vector", {
  a <- c()
  expect_that(is.empty(a), equals(TRUE))
})

test_that("non-empty list", {
  a <- list(a=1,2,3)
  expect_that(is.empty(a), equals(FALSE))
})
test_that("empty list", {
  a <- list()
  expect_that(is.empty(a), equals(TRUE))
})

test_that("non-empty data.frame", {
  a <- data.frame(a=1:3,b=2,c=3)
  expect_that(is.empty(a), equals(FALSE))
})
test_that("empty data.frame", {
  a <- data.frame(a=NULL, b=NULL)
  expect_that(is.empty(a), equals(TRUE))
})

