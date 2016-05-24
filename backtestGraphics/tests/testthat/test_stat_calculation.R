context("Test for the function 'stat_calculation'")

load("test_stat_calculation.RData")

test_that("stat_calculation function", {
  result.1 <- backtestGraphics:::stat_calculation(x.list)
  
  expect_equal(result.1, truth.1, label = "Failed the test for calculating summary statistics")
})