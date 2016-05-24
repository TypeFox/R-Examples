library(functools)
context("Reject()")

less_than_5 <- function(x) x < 5
more_than_m <- function(x) x > "m"
foo <- 1:10
test_that("Produces the correct output.", {
  expect_equal(Reject(less_than_5, foo), Filter(Negate(less_than_5), foo))
  expect_equal(Reject(more_than_m, letters), Filter(Negate(more_than_m), letters))
})

test_that("Produces the correct output type.", {
  expect_is(Reject(less_than_5, foo), "integer")
  expect_is(Reject(more_than_m, letters), "character")
})

test_that("Produces the correct errors.", {
  expect_error(Reject(foo, less_than_5))
  expect_error(Reject(letters, more_than_m))
})
