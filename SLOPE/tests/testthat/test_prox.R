test_that("prox_sorted_L1 agrees with isotone package", {
  if (!requireNamespace('isotone', quietly=TRUE))
    skip('isotone package not available')
  n = 20
  mu = 1.5 * (n:1)
  y = sort(abs(rnorm(n, mean=mu)), decreasing=TRUE)
  lambda = n:1
  expect_equal(prox_sorted_L1_C(y,lambda), prox_sorted_L1_isotone(y,lambda))
})