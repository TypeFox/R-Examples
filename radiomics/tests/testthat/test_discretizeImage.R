library(radiomics)
context("discretizeImage")

test_that("Image discretization ", {
  expect_error(discretizeImage(data.frame()))
  expect_error(discretizeImage())
  expect_error(discretizeImage(c()))
  expect_error(discretizeImage(list()))
  expect_warning(discretizeImage(matrix()))
  
})
