context("GetSizeOfGDELT")

test_that("returns valid size for GDELT", {
  expect_true(GetSizeOfGDELT() > 8.94)  # as of 1/16/2014
})
