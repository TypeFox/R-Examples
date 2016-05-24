library(generator)
library(detector)
context("email addresses")

test_that("Produces the correct output.", {
  expect_equal(is_email_address("hello"), FALSE)
  expect_equal(is_email_address("hello@world.edu"), TRUE)
  expect_equal(has_email_addresses("hello"), FALSE)
  expect_equal(has_email_addresses("hello@world.edu"), TRUE)
})

test_that("Produces the correct output type.", {
  expect_is(is_email_address("hello"), "logical")
  expect_is(is_email_address("hello@world.edu"), "logical")
  expect_is(has_email_addresses("hello"), "logical")
  expect_is(has_email_addresses("hello@world.edu"), "logical")
})

test_that("Produces the correct errors.", {
})

