context("linestring")

test_that("convert linestring works", {
  # empty
  empty <- linestring("empty")
  expect_is(empty, "character")
  expect_equal(empty, "LINESTRING EMPTY")

  # numeric
  b <- linestring(c(100.000, 3.101), c(101.000, 2.100), c(3.140, 2.180), fmt = 0)
  expect_is(b, "character")
  expect_match(b, "LINESTRING")
  expect_equal(b, "LINESTRING (100.000 3.101, 101.0 2.1, 3.14 2.18)")

  # data.frame
  str <- "LINESTRING (-116.4 45.2, -118.0 47.0)"
  df <- data.frame(lon = c(-116.4,-118), lat = c(45.2,47))
  lsdf <- linestring(df, fmt = 1)
  expect_is(lsdf, "character")
  expect_match(lsdf, "LINESTRING")
  expect_equal(lsdf, str)

  # matrix
  mat <- matrix(c(df$lon, df$lat), ncol = 2)
  lsmat <- linestring(mat, fmt = 1)
  expect_is(lsmat, "character")
  expect_match(lsmat, "LINESTRING")
  expect_equal(lsmat, str)

  # list
  lslist <- linestring(list(c(100.000, 0.000), c(101.000, 1.000)), fmt = 0)
  expect_is(lslist, "character")
  expect_match(lslist, "LINESTRING")
  expect_equal(lslist, "LINESTRING (100 0, 101 1)")
})

test_that("linestring fails correctly", {
  expect_error(linestring(-116.4), "LINESTRING input should be of length 2")
  expect_error(linestring(), "no applicable method")
  expect_error(linestring(NA), "no applicable method")
  expect_error(linestring("a", "Adf", mtcars), "The following strings are not WKT")
})
