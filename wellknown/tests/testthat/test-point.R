context("point")

test_that("convert point works", {
  # empty point
  empty <- point("empty")
  expect_is(empty, "character")
  expect_equal(empty, "POINT EMPTY")

  ## numeric, single point
  pt_a <- point(-116.4, 45.2)
  pt_b <- point(0, 1)
  expect_is(pt_a, "character")
  expect_is(pt_b, "character")
  expect_match(pt_a, "POINT")
  expect_match(pt_b, "POINT")
  expect_equal(point(-116.4, 45.2), "POINT (-116.4000000000000057 45.2000000000000028)")
  expect_equal(point(-116.4, 45.2, fmt = 1), "POINT (-116.4 45.2)")

  ## single point, from data.frame
  df <- data.frame(lon = -116.4, lat = 45.2)
  df_a <- point(df, fmt = 2)
  expect_is(df_a, "character")
  expect_is(df_a[1], "character")
  expect_equal(df_a[1], "POINT (-116.40 45.20)")

  ## many points, from data.frame
  df2 <- us_cities[1:5,c('lat','long')]
  df_b <- point(df2)
  expect_is(df_b, "character")
  expect_is(df_b[1], "character")
  expect_equal(df_b[1], "POINT (32.4500000000000028 -99.7399999999999949)")

  ## matrix
  ussmall <- us_cities[1:5, ]
  df <- data.frame(long = ussmall$long, lat = ussmall$lat)
  mat <- matrix(c(df$long, df$lat), ncol = 2)
  ptmat <- point(mat, fmt = 0)
  expect_is(ptmat[1], "character")
  expect_is(ptmat[2], "character")
  expect_equal(ptmat[1], "POINT (-99.74 32.45)")

  ## single point, from a list
  ls_a <- point(list(c(100.0, 3.1)), fmt = 2)
  expect_is(ls_a, "character")
  expect_is(ls_a[1], "character")
  expect_equal(ls_a[1], "POINT (100.00 3.10)")
})

test_that("point fails correctly", {
  expect_error(point(-116.4), "POINT input should be of length 2")
  expect_error(point(), "no applicable method")
  expect_error(point(NA), "no applicable method")
  expect_error(point("a", "Adf"), "The following strings are not WKT")
})
