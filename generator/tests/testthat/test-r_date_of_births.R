library(generator)
context("r_date_of_births()")

test_that("Produces the correct output.", {
  expect_equal(length(r_date_of_births(100)), 100)
})

test_that("Produces the correct output type.", {
  expect_is(r_date_of_births(100), "Date")
})

test_that("Produces the correct errors.", {
})

