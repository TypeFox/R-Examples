test_that("AD Test Stat.",{

  # Success cases were/can be tested by comparing plots and confidence intervals.
  
  # Error - 1
  expect_error(val <- ADTestStat(1000, 100, 1.05))

  # Error - 2
  expect_error(val <- ADTestStat(1000, 100, -.95))
  
})