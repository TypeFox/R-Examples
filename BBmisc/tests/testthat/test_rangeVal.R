
context("rangeVal")

test_that("rangeVal", {
  expect_equal(rangeVal(c(1, 5)), 4)
  expect_equal(rangeVal(1), 0)
  expect_equal(rangeVal(1:3), 2)
  
  # NAs
  expect_equal(rangeVal(c(1, 2, NA)), NA_real_)
  expect_equal(rangeVal(c(1, 2, NA), na.rm = TRUE), 1)
  expect_equal(rangeVal(c(NA_real_), na.rm = TRUE), NA_real_)
})
