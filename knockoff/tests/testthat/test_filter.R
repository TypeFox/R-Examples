test_that('knockoff.filter verifies input dimensions', {
  expect_error(knockoff.filter(rnorm_matrix(10, 10), rnorm(10)), 'dimensions')
  expect_error(knockoff.filter(rnorm_matrix(20, 10), rnorm(19)))
  expect_warning(knockoff.filter(rnorm_matrix(20, 15), rnorm(20)), 'dimensions')
})