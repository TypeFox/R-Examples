library(generator)
context("r_full_names()")

test_that("Produces the correct output.", {
  expect_equal(length(r_full_names(100)), 100)
})

test_that("Produces the correct output type.", {
  expect_is(r_full_names(100), "character")
})

test_that("Produces the correct errors.", {
})

