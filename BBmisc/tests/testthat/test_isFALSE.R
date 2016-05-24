context("isFALSE")

test_that("isFALSE", {
  expect_equal(isFALSE(FALSE), TRUE)
  expect_equal(isFALSE(TRUE), FALSE)
  expect_equal(isFALSE(0), FALSE)
})  
