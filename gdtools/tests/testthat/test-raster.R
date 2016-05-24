context("raster")

test_that("raster_str", {
  r <- as.raster(matrix(hcl(0, 80, seq(40, 80, 10)), nrow = 5, ncol = 4))
  code <- raster_str(r, width = 50, height = 50)
  expect_more_than(nchar(x = code), 50)
})

test_that("raster_write", {
  r <- as.raster(matrix(hcl(0, 80, seq(40, 80, 10)), nrow = 5, ncol = 4))
  filename <- tempfile(fileext = ".png")
  raster_write(r, path = filename, width = 50, height = 50)
  expect_true(file.exists(filename))
})
