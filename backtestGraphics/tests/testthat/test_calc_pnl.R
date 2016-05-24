context("Test for function 'calc_pnl'")

load("test_calc_pnl.RData")

test_that("calc_pnl function", {
  
  result.1 <- backtestGraphics:::calc_pnl(x)
  
  expect_equal(result.1, truth.1, label = "Failed the test to generate information about returns and profits")
  
})