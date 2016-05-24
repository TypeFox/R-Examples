context("Test exif")

test_that("Error handling works", {
  error_file <- system.file("extdata/shark_cat_test_img.jpg", package="exif")
  testthat::expect_error(read_exif(error_file))
})

test_that("Single files can be read", {
  single_file <- system.file("extdata/dog_test_img.jpg", package="exif")
  results <- read_exif(single_file)
  expect_equal(nrow(results),1)
  expect_equal(ncol(results),31)
})