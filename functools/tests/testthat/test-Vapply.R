library(functools)
context("Vapply()")

x <- runif(10000)
test_that("Produces the correct output.", {
  expect_equal(Vapply(x, function(x) x * x, numeric(1)),
               vapply(x, function(x) x * x, numeric(1)))
})

test_that("Produces the correct output type.", {
  expect_is(Vapply(x, Identity, numeric(1)), "numeric")
})

test_that("Produces the correct errors.", {
  expect_error(Vapply(Identity, x))
})
