library(functools)
context("Apply()")

x <- matrix(runif(10000), ncol = 5)
test_that("Produces the correct output.", {
  expect_equal(Apply(x, mean, 1), apply(x, 1, mean))
  expect_equal(Apply(x, mean, 2), apply(x, 2, mean))
  expect_equal(Apply(x, sum, 1), apply(x, 1, sum))
  expect_equal(Apply(x, sum, 2), apply(x, 2, sum))
})

test_that("Produces the correct output type.", {
  expect_is(Apply(x, mean, 1), "numeric")
  expect_is(Apply(x, mean, 2), "numeric")
  expect_is(Apply(x, sum, 1), "numeric")
  expect_is(Apply(x, sum, 2), "numeric")

})

test_that("Produces the correct errors.", {
  expect_error(Apply(x, 1, mean))
  expect_error(Apply(x, 2, mean))
})
