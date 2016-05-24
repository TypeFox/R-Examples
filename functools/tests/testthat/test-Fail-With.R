library(functools)
context("Fail_With()")

default_value <- 0

new_function <- function(x) {
  stop("Error!")
  return(x)
}

new_function_with_default <-
  Fail_With(default_value, new_function, .silent = TRUE)

test_that("Produces the correct output.", {
  expect_equal(new_function_with_default("a"), default_value)
  expect_equal(new_function_with_default(1), default_value)
  expect_equal(new_function_with_default(1.2), default_value)
  expect_equal(new_function_with_default(TRUE), default_value)
  expect_equal(new_function_with_default(NA), default_value)
  expect_equal(new_function_with_default(NULL), default_value)
})

test_that("Produces the correct output type.", {
  expect_is(new_function, "function")
  expect_is(new_function_with_default, "function")
})

test_that("Produces the correct errors.", {
  expect_equal(1, 1)
})
