library(functools)
context("Partial()")

new_mean <- Partial(mean, na.rm = TRUE)

x <- c(NA, 1:10)

test_that("Produces the correct output.", {
  expect_equal(new_mean(x), mean(x, na.rm = TRUE))
})

test_that("Produces the correct output type.", {
  expect_is(new_mean, "function")
})

test_that("Produces the correct errors.", {
  expect_error(Partial()())
})
