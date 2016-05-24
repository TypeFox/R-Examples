test_that("Binomial Backtest Works.",{
  # Success - 1
  # x < n * p
  expect_equal(0.7358, BinomialBacktest(1, 100, 0.99), tolerance=0.001)
  
  # Success - 2
  # x >= n * p
  expect_equal(0, BinomialBacktest(91, 100, 0.90), tolerance=0.001)
  
  # Error - 1
  expect_error(val <- BinomialBacktest(35, 30, 0.95))
  
  # Error - 2
  expect_error(val <- BinomialBacktest(5, 30, 1.5))
  
  # Error - 3
  expect_error(val <- BinomialBacktest(5, 30, -.95))

})