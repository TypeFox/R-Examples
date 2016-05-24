context("wktview")

skip_if_not_installed <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE))
    return()
  testthat::skip(paste0(pkg, " not installed"))
}

# from point  ----------------
test_that("wktview works with point", {
  skip_if_not_installed("leaflet")

  str <- "POINT (-116.4000000000000057 45.2000000000000028)"
  pt1 <- wktview(str)
  pt2 <- wktview(point(list(c(100.0, 3.101))))

  expect_is(pt1, c("leaflet", "htmlwidget"))
  expect_is(pt1$x, "list")
  expect_null(pt1$x$fitBounds)
  expect_is(pt1$x$calls, "list")
  expect_is(pt2, c("leaflet", "htmlwidget"))
})

# from multipoint -----------------------------
test_that("wktview works with multipoint", {
  skip_if_not_installed("leaflet")

  mpt <- multipoint(c(100.000, 3.101), c(101.000, 2.100), c(3.140, 2.180))
  mpt1 <- wktview(mpt)

  expect_is(mpt1, c("leaflet", "htmlwidget"))
  expect_is(mpt1$x, "list")
  expect_is(mpt1$x$calls, "list")
})

# from polygon ----------------
test_that("wktview works with polygon input", {
  skip_if_not_installed("leaflet")

  ply <- polygon(c(100.001, 0.001), c(101.12345, 0.001), c(101.001, 1.001), c(100.001, 0.001), fmt = 2)
  poly1 <- wktview(ply)

  expect_is(poly1, c("leaflet", "htmlwidget"))
  expect_is(poly1$x, "list")
  expect_null(poly1$x$fitBounds)
  expect_is(poly1$x$calls, "list")
})

# from multipolygon  --------------------------------
test_that("wktview works with multipolygon", {
  skip_if_not_installed("leaflet")

  df <- data.frame(long = c(30, 45, 10, 30), lat = c(20, 40, 40, 20))
  df2 <- data.frame(long = c(15, 40, 10, 5, 15), lat = c(5, 10, 20, 10, 5))
  mpoly <- multipolygon(df, df2, fmt = 0)
  mpoly1 <- wktview(mpoly)

  expect_is(mpoly1, c("leaflet", "htmlwidget"))
  expect_is(mpoly1$x, "list")
  expect_is(mpoly1$x$calls, "list")
})

# from linestring  ----------------
test_that("wktview works with linestring", {
  skip_if_not_installed("leaflet")

  line <- linestring("LINESTRING (-116.4 45.2, -118.0 47.0)")
  line1 <- wktview(line)

  expect_is(line1, c("leaflet", "htmlwidget"))
  expect_is(line1$x, "list")
  expect_null(line1$x$fitBounds)
  expect_is(line1$x$calls, "list")
})

# from multilinestring  ----------------
test_that("wktview works with multilinestring", {
  skip_if_not_installed("leaflet")

  df <- data.frame(long = c(30, 45, 10), lat = c(20, 40, 40))
  df2 <- data.frame(long = c(15, 40, 10), lat = c(5, 10, 20))
  mline <- multilinestring(df, df2, fmt=0)
  mline1 <- wktview(mline)

  expect_is(mline1, c("leaflet", "htmlwidget"))
  expect_is(mline1$x, "list")
  expect_null(mline1$x$fitBounds)
  expect_is(mline1$x$calls, "list")
})

# from geometrycollection  ----------------
# leaflet is not working with geometrycollection yet
# seems like it should though cause leaflet supports the geometry
test_that("wktview works with geometrycollection", {
#   skip_if_not_installed("leaflet")
#
#   gc <- geometrycollection(
#     point(-116.4, 45.2),
#     linestring("LINESTRING (-116.4 45.2, -118.0 47.0)")
#   )
#
#   expect_is(wktview(gc), c("leaflet", "htmlwidget"))
})

# wktview fails well  ----------------
test_that("wktview fails well", {
  skip_if_not_installed("leaflet")

  # missing arguments
  expect_error(wktview(), "no applicable method")
  # wrong input
  expect_error(wktview("a"))
  # no method for data.frame's
  expect_error(wktview(mtcars), "no applicable method")
})
