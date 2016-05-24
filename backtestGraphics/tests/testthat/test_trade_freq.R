context("Test for the function 'trade_freq'")

load("test_trade_freq.RData")

test_that("trade_ferq function", {
  result.1 <- backtestGraphics:::trade_freq(x)
  expect_equal(truth.1, result.1, label = "Failed the test for trading frequency")
})