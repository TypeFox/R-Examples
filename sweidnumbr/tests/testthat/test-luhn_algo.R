
context("luhn_algo")

test_that("fails", {
  expect_that(luhn_algo(123), throws_error())
  expect_that(luhn_algo("850504333", list(1)), throws_error())
  expect_that(luhn_algo(123456789, c(2, 1, 2, 1, 2, 1, 2, 1, 2, 0)), throws_error())
  expect_that(luhn_algo(c(198505043334, 1212121212)), throws_error())
})

test_that("1 pin works", {
  expect_that(suppressMessages(luhn_algo(198505043334)), is_equivalent_to(4))
  expect_that(suppressMessages(luhn_algo(19850504333)), is_equivalent_to(4))
  expect_that(suppressMessages(luhn_algo(8505043334)), is_equivalent_to(4))
  expect_that(suppressMessages(luhn_algo(850504333)), is_equivalent_to(4))
  expect_that(suppressMessages(luhn_algo("8505043334")), is_equivalent_to(4))
  expect_that(luhn_algo(8505043334), shows_message())
  expect_that(suppressMessages(luhn_algo("8505043334", multiplier = c(2, 1, 2, 1, 2, 1, 2, 1, 2, 0))), 
              is_equivalent_to(4))  
})

test_that(desc="Handle NA in luhn_algo",{
  expect_true(is.na(luhn_algo(id = c(NA,198501169885))[1]))
})

test_that("multiple pin works", {
  expect_that(suppressMessages(luhn_algo(c(198505043334,121212121212))), is_equivalent_to(c(4, 2)))
})