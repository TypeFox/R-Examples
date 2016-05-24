context("polygon")

test_that("polygon works", {
  ## empty
  empty <- polygon("empty")
  expect_is(empty, "character")
  expect_equal(empty, "POLYGON EMPTY")

  ## numeric
  poly <- polygon(c(100, 0.1), c(101, 0.1), c(101, 1.1), c(100, 0.1), fmt = 0)
  expect_is(poly, "character")
  expect_match(poly, "POLYGON")
  expect_equal(poly, "POLYGON ((100.0 0.1, 101.0 0.1, 101.0 1.1, 100.0 0.1))")

  ## from data.frame
  str <- "POLYGON ((-81.52 41.08, -122.26 37.77, -84.18 31.58, -73.80 42.67, -81.52 41.08))"
  df <- us_cities[2:5,c('long','lat')]
  df <- rbind(df, df[1, ])
  polydf <- polygon(df, fmt = 2)
  expect_is(polydf, "character")
  expect_equal(polydf, str)

  ## from data.frame, many
  str <- "POLYGON ((-81.52 41.08, -122.26 37.77, -84.18 31.58, -73.80 42.67, -81.52 41.08), (-85.9 37.5, -85.9 35.3, -93.0 35.3, -93.0 37.5, -85.9 37.5))"
  df2 <- data.frame(long = c(-85.9, -85.9, -93, -93, -85.9),
                    lat = c(37.5, 35.3, 35.3, 37.5, 37.5))
  manypolydf <- polygon(df, df2, fmt = 0)
  expect_is(manypolydf, "character")
  expect_equal(manypolydf, str)

  ## single polygon, from a matrix
  str <- "POLYGON ((-81.52 41.08, -122.26 37.77, -84.18 31.58, -73.80 42.67, -81.52 41.08))"
  mat <- matrix(c(df$long, df$lat), ncol = 2)
  polymat <- polygon(mat, fmt = 0)
  expect_is(polymat, "character")
  expect_equal(polymat, str)

  ## from a list
  str <- "POLYGON ((100.001 0.001, 101.1235 0.0010, 101.001 1.001, 100.001 0.001))"
  ply <- list(c(100.001, 0.001), c(101.12345, 0.001), c(101.001, 1.001), c(100.001, 0.001))
  polylist <- polygon(ply, fmt = 2)
  expect_is(polylist, "character")
  expect_equal(polylist, str)
})

test_that("polygon fails correctly", {
  expect_error(polygon(-116.4), "POLYGON input should be of length 2")
  expect_error(polygon(), "no applicable method")
  expect_error(polygon(NA), "no applicable method")
  expect_error(polygon("a", "Adf"), "The following strings are not WKT")
})
