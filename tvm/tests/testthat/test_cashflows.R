context("Cashflows")

test_that("Cashflow for a bullet loan", {
  l <- loan(rate = 0.1, maturity = 4, amt = 1, type = "bullet", grace_int = 0, grace_amort = 0)
  expect_that(l$cf, equals(c(0.1, 0.1, 0.1, 1.1)))
})

test_that("Cashflow for a german loan", {
  l <- loan(rate = 0.1, maturity = 4, amt = 1, type = "german", grace_int = 0, grace_amort = 0)
  expect_that(l$cf, equals(0.25 + c(1, 0.75, 0.5, 0.25) * 0.1))
})

test_that("Cashflow for a french loan", {
  l <- loan(rate = 0.1, maturity = 4, amt = 1, type = "french", grace_int = 0, grace_amort = 0)
  pmt <- 1 / sum (1 / (rep_len(1 + 0.1, 4) ^ (1:4)))
  expect_that(l$cf, equals(rep_len(pmt, 4)))
})