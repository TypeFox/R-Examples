context("multipoint")

test_that("convert multipoint works", {
  # empty multipoint
  empty <- multipoint("empty")
  expect_is(empty, "character")
  expect_equal(empty, "MULTIPOINT EMPTY")

  # numeric
  b <- multipoint(c(100.000, 3.101), c(101.000, 2.100), c(3.140, 2.180), fmt = 0)
  expect_is(b, "character")
  expect_match(b, "MULTIPOINT")
  expect_equal(b, "MULTIPOINT ((100.000 3.101), (101.0 2.1), (3.14 2.18))")

  # data.frame
  str <- "MULTIPOINT ((-99.74 32.45), (-81.52 41.08), (-122.26 37.77), (-84.18 31.58), (-73.80 42.67))"
  df <- us_cities[1:5, c('long', 'lat')]
  multidf <- multipoint(df, fmt = 0)
  expect_is(multidf, "character")
  expect_match(multidf, "MULTIPOINT")
  expect_equal(multidf, str)

  # matrix
  mat <- matrix(c(df$long, df$lat), ncol = 2)
  multimat <- multipoint(mat, fmt = 0)
  expect_is(multimat, "character")
  expect_match(multimat, "MULTIPOINT")
  expect_equal(multimat, str)

  # list
  multilist <- multipoint(list(c(100.000, 3.101), c(101.000, 2.100), c(3.140, 2.180)), fmt = 0)
  expect_is(multilist, "character")
  expect_match(multilist, "MULTIPOINT")
  expect_equal(multilist, "MULTIPOINT ((100.000 3.101), (101.0 2.1), (3.14 2.18))")
})

test_that("multipoint fails correctly", {
  expect_error(multipoint(-116.4), "POINT input should be of length 2")
  expect_error(multipoint(), "no applicable method")
  expect_error(multipoint(NA), "no applicable method")
  expect_error(multipoint("a", "Adf"), "The following strings are not WKT")
})
