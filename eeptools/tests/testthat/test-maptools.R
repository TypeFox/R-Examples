# Map functions
context("Does mapmerge work")

test_that("Map merge matches line up correctly", {
  xx <- maptools::readShapePoly(system.file("shapes/sids.shp", package="maptools")[1], 
                                IDvar="FIPSNO")
  yy <- as(xx, "data.frame")
  yy$newvar <- sample(letters, nrow(yy), replace=TRUE)
  yy <- subset(yy, select=c("FIPS", "newvar"))
  newpoly <- mapmerge(xx, yy, xid="FIPS", yid="FIPS")
  expect_identical(newpoly@data[, c("FIPS", "newvar")], yy)
  expect_identical(nrow(xx@data), nrow(newpoly@data))
})

test_that("Map merge results in an object with only data modified", {
  xx <- maptools::readShapePoly(system.file("shapes/sids.shp", package="maptools")[1], 
                                IDvar="FIPSNO")
  yy <- as(xx, "data.frame")
  yy$newvar <- sample(letters, nrow(yy), replace=TRUE)
  yy <- subset(yy, select=c("FIPS", "newvar"))
  newpoly <- mapmerge(xx, yy, xid="FIPS", yid="FIPS")
  expect_is(newpoly, "SpatialPolygonsDataFrame")
  expect_false(identical(xx@data, newpoly@data))
  expect_identical(xx@polygons, newpoly@polygons)
  expect_identical(xx@proj4string, newpoly@proj4string)
  expect_identical(xx@plotOrder, newpoly@plotOrder)
  expect_identical(xx@bbox, newpoly@bbox)
})

context("Test that ggmapmerge works")

test_that("It returns the correct data.frame", {
  xx <- maptools::readShapePoly(system.file("shapes/sids.shp", package="maptools")[1], IDvar="FIPSNO")
  plotobj <- ggmapmerge(xx, "FIPS")
  expect_is(plotobj, "data.frame")
  expect_identical(names(plotobj), c("id", "long", "lat", "order", "hole", 
                                     "piece", "group", "AREA", "PERIMETER", 
                                     "CNTY_", "CNTY_ID", "NAME", "FIPSNO", 
                                     "CRESS_ID", "BIR74", "SID74", "NWBIR74", 
                                     "BIR79", "SID79", "NWBIR79"))
  expect_equal(nrow(plotobj), 2529)
})

