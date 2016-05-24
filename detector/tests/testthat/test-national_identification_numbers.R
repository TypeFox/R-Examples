library(generator)
library(detector)
context("national identification numbers")

test_that("Produces the correct output.", {
  expect_equal(is_national_identification_number("hello"), FALSE)
  expect_equal(is_national_identification_number("000-00-0000"), TRUE)
  expect_equal(has_national_identification_numbers("hello"), FALSE)
  expect_equal(has_national_identification_numbers("000-00-0000"), TRUE)
})

test_that("Produces the correct output type.", {
  expect_is(is_national_identification_number("hello"), "logical")
  expect_is(is_national_identification_number("000-00-0000"), "logical")
  expect_is(has_national_identification_numbers("hello"), "logical")
  expect_is(has_national_identification_numbers("000-00-0000"), "logical")
})

test_that("Produces the correct errors.", {
})

