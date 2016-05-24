library(functools)
context("Tapply()")

groups <- as.factor(rbinom(size = 32, n = 10000, prob = 0.1))
n <- 17; fac <- factor(rep(1:3, length = n), levels = 1:5)

test_that("Produces the correct output.", {
  expect_equal(Tapply(groups, groups, .f = length), tapply(groups, groups, length))
})

test_that("Produces the correct output type.", {
  expect_is(Tapply(groups, groups, .f = length), "array")
})

test_that("Produces the correct errors.", {
  expect_error(Tapply(length, groups, groups))
})
