context("Test getCountries")

test_that("Testing functionality", {
  skip_on_cran()
  expect_is(test <- getCountries(), 'data.frame')
  expect_equal(ncol(test),2)
})