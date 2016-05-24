library(generator)
context("r_latitudes()")

test_that("Produces the correct output.", {
  expect_equal(length(r_latitudes(100)), 100)
})

test_that("Produces the correct output type.", {
  expect_is(r_latitudes(100), "numeric")
})

test_that("Produces the correct errors.", {
})

