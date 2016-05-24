library(glmmsr)
context("Utility functions")

test_that("correctly detects whether a formula has random effects", {
  expect_equal(has_re(y ~ 0 + (1 | p)), TRUE)
  expect_equal(has_re(y ~ x), FALSE)
  expect_equal(has_re(y ~ x*a + b), FALSE)
  expect_equal(has_re(y ~ (1 | p) + x), TRUE)
  expect_equal(has_re(y ~ (1 | p) + (1 | q)), TRUE)
  expect_equal(has_re(y ~ (p | q)), TRUE)
  expect_equal(has_re(y ~ (p || q)), TRUE)
})
