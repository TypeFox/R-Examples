library(functools)
context("Sapply()")

x <- list(a = runif(10000),
          beta = rexp(10000),
          logic = sample(x = c(TRUE, FALSE), size = 10000, replace = TRUE))
test_that("Produces the correct output.", {
  expect_equal(Sapply(x, quantile), sapply(x, quantile))
  expect_equal(Sapply(x, min), sapply(x, min))
  expect_equal(Sapply(x, max), sapply(x, max))
})

test_that("Produces the correct output type.", {
  expect_is(Sapply(x, quantile), "matrix")
  expect_is(Sapply(x, min), "numeric")
  expect_is(Sapply(x, max), "numeric")
})

test_that("Produces the correct errors.", {
  expect_error(Sapply(quantile, x))
  expect_error(Sapply(min, x))
  expect_error(Sapply(max, x))
})
