context("chop")

suppressPackageStartupMessages(library("rgeos"))
wkt <- "POLYGON((-180 -20, -140 55, 10 0, -140 -60, -180 -20))"
poly <- readWKT(wkt)

test_that("chop inputs", {
  expect_is(wkt, "character")
  expect_match(wkt, "POLYGON")
  expect_is(poly, "SpatialPolygons")
})

test_that("SpatialPolygons input works", {
  aa <- chop(poly)

  expect_is(aa, "SpatialPolygons")
  expect_is(aa@polygons[[1]], "Polygons")
  expect_is(aa@proj4string, "CRS")
  expect_equal(length(aa@polygons), 139)
})

test_that("WKT input works", {
  aa <- chop(wkt)

  expect_is(aa, "SpatialPolygons")
  expect_is(aa@polygons[[1]], "Polygons")
  expect_is(aa@proj4string, "CRS")
  expect_equal(length(aa@polygons), 139)
})

test_that("geojson character input works", {
  file <- system.file("examples", "sample1.geojson", package = "geoaxe")
  x <- readLines(file)
  aa <- chop(x)

  expect_is(x, "character")
  expect_match(x, "coordinates")
  expect_is(aa, "SpatialPolygons")
  expect_is(aa@polygons[[1]], "Polygons")
  expect_is(aa@proj4string, "CRS")
  expect_equal(length(aa@polygons[[1]]), 1)
})
