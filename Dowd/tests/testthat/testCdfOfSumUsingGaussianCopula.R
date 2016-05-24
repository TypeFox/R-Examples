test_that("Cdf Of Sum Using Gaussian Copula.",{
  # Success - 1
  expect_equal(0.3446, CdfOfSumUsingGaussianCopula(0.9, 0.3, 1.5, 1.2, 1.5, 0.6, 25), tolerance=0.1)
  
  # Success - 2
  expect_equal(0.9939, CdfOfSumUsingGaussianCopula(0.1, -0.3, -3, .2, 1.5, 0.6, 10), tolerance=0.1)
  
})
