library(generator)
context("r_national_identification_numbers()")

test_that("Produces the correct output.", {
  expect_equal(length(r_national_identification_numbers(100)), 100)
})

test_that("Produces the correct output type.", {
  expect_is(r_national_identification_numbers(100), "character")
})

test_that("Produces the correct errors.", {
})

