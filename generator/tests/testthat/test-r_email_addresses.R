library(generator)
context("r_email_addresses()")

test_that("Produces the correct output.", {
  expect_equal(length(r_email_addresses(100)), 100)
})

test_that("Produces the correct output type.", {
  expect_is(r_email_addresses(100), "character")
})

test_that("Produces the correct errors.", {
})

