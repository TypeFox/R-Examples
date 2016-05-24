context("as_SpatialPolgyons")

file1 <- system.file("examples", "sample1.geojson", package = "geoaxe")
file2 <- system.file("examples", "sample2.geojson", package = "geoaxe")
file3 <- system.file("examples", "sample3.geojson", package = "geoaxe")

wo_feat <- jsonlite::fromJSON(file1, FALSE)
w_feat <- jsonlite::fromJSON(file2, FALSE)
many_feat <- jsonlite::fromJSON(file3, FALSE)

wo_feat_json <- jsonlite::toJSON(wo_feat, auto_unbox = TRUE)
w_feat_json <- jsonlite::toJSON(w_feat, auto_unbox = TRUE)
many_feat_json <- jsonlite::toJSON(many_feat, auto_unbox = TRUE)

test_that("as_SpatialPolgyons works for list input", {
  aa <- as_SpatialPolygons(wo_feat)
  bb <- as_SpatialPolygons(w_feat)
  cc <- suppressMessages(as_SpatialPolygons(many_feat))

  expect_is(aa, "SpatialPolygons")
  expect_is(aa@polygons[[1]], "Polygons")
  expect_is(aa@proj4string, "CRS")
  expect_equal(length(aa@polygons[[1]]@Polygons), 1)

  expect_is(bb, "SpatialPolygons")
  expect_is(bb@polygons[[1]], "Polygons")
  expect_is(bb@proj4string, "CRS")
  expect_equal(length(bb@polygons[[1]]@Polygons), 1)

  expect_is(cc, "SpatialPolygons")
  expect_is(cc@polygons[[1]], "Polygons")
  expect_is(cc@proj4string, "CRS")
  expect_equal(length(cc@polygons[[1]]@Polygons), 2)
  expect_message(as_SpatialPolygons(many_feat),
                 "only polygons supported")
})
