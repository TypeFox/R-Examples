context("Test for function 'best_worst_three'")

load("test_best_worst_three.RData")

test_that("best_worst_three function", {
  
  result.1 <- backtestGraphics:::best_worst_three(x)
  
  expect_equal(result.1, truth.1, label = "Failed the test to find the best/worst three performers")

})