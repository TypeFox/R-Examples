test_that("Cdf Of Sum Using Product Copula.",{
  # Success - 1
  expect_equal(0.1628, CdfOfSumUsingProductCopula(0.9, .3, 2.1, .3, 1.5), tolerance=0.001)
  
  # Success - 2
  expect_equal(0.5839, CdfOfSumUsingProductCopula(0.1, 1.5, -2.5, 0.4, 5.1), tolerance=0.001)
  
})