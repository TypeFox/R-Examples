context("Currency conversion")

test_that("Simple amounts can be converted",{
  expect_that(to_currency(12), equals("£12.00"))
})

test_that("Complex amounts can be converted",{
  expect_that(to_currency(12000), equals("£12,000.00"))
})

test_that("Complex amounts can be converted with non-default arguments",{
  expect_that(to_currency(12000, currency_symbol = "$", symbol_first = FALSE,
                          decimal_size = 3, decimal_delim = "q"), equals("12,000q000$"))
})

test_that("Simple strings can be converted into amounts", {
  expect_that(from_currency("$120,000", "."), equals(120000.00))
})

test_that("Complex strings can be converted into amounts", {
  expect_that(from_currency("$120,000.022", "."), equals(120000.022))
})
