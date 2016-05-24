test_that("JarqueBera Backtest.",{
  # Success
  expect_equal(0.6940, JarqueBeraBacktest(-0.0758, 2.8888, 500), tolerance=0.01)
  expect_equal(1, JarqueBeraBacktest(0, 3, 100), tolerance=0.01)
  
})