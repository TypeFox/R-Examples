test_that("SLOPE solver agrees with TFOCS solver in MATLAB", {
  skip_on_cran()
  if (!requireNamespace('R.matlab', quietly=TRUE))
    skip('R.matlab package not available')

  problem = random_problem(100, 50, sigma=1)
  X = problem$X
  y = problem$y
  result.R = SLOPE(X, y, sigma=1, solver='default')
  result.matlab = SLOPE(X, y, sigma=1, solver='matlab')

  expect_equal(result.R$selected, result.matlab$selected)
  expect_equal(result.R$beta, result.matlab$beta, tol=1e-6)
})