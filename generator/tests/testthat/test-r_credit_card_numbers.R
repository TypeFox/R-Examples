library(generator)
context("r_credit_card_numbers()")

test_that("Produces the correct output.", {
  expect_equal(length(r_credit_card_numbers(100)), 100)
})

test_that("Produces the correct output type.", {
  expect_is(r_credit_card_numbers(100), "character")
})

test_that("Produces the correct errors.", {
})

