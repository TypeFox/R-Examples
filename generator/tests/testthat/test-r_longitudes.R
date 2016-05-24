library(generator)
context("r_longitudes()")

test_that("Produces the correct output.", {
  expect_equal(length(r_longitudes(100)), 100)
})

test_that("Produces the correct output type.", {
  expect_is(r_longitudes(100), "numeric")
})

test_that("Produces the correct errors.", {
})

