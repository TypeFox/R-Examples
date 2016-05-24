test_that("Cdf Of Sum Using Gumbel Copula.",{
  # Success - 1
  expect_equal(0.0259, CdfOfSumUsingGumbelCopula(0.1, 1.3, 2.5, 1.2, 1.5, 1.2), tolerance=0.1)
  
  # Success - 2
  expect_equal(0.8896, CdfOfSumUsingGumbelCopula(0.1, 1.5, -2.5, 0.4, 5.1, 2.1), tolerance=0.1)
  
  # Error - 1
  expect_error(val <- CdfOfSumUsingGumbelCopula(0.1, 1.5, -2.5, 0.4, 5.1, .9))
})
