context("Test for function 'best_worst_month'")

load("test_best_worst_month.RData")

test_that("best_worst_month function", {
  
  result.1 <- backtestGraphics:::best_worst_month(x)
  
  expect_equal(result.1, truth.1, label = "Failed the test to calculate the best/worst months")
})
