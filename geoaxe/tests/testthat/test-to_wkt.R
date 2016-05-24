context("to_wkt")

suppressPackageStartupMessages(library("rgeos"))
wkt <- "POLYGON((-180 -20, -140 55, 10 0, -140 -60, -180 -20))"
poly <- readWKT(wkt)
polys <- chop(x = poly)

test_that("to_wkt works for SpatialPolygons input", {
  aa <- to_wkt(polys)

  expect_is(aa, "axewkt")
  expect_is(aa[[1]], "character")
  expect_match(aa[[1]], "POLYGON")
  expect_equal(length(aa), 139)
})

test_that("to_wkt works for axewkt (self) input", {
  aa <- to_wkt(polys)
  bb <- to_wkt(aa)

  expect_is(aa, "axewkt")
  expect_is(bb, "axewkt")
  expect_identical(aa, bb)
})
