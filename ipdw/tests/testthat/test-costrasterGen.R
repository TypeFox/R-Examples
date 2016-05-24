context("costrasterGen")

test_that("returns correct class and warn on bad projection", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()
  
  Sr1 <- Polygon(cbind(c(0, 0, 1, 1, 0), c(0, 12, 12, 0, 0)))
  Sr4 <- Polygon(cbind(c(9, 9, 10, 10, 9), c(0, 12, 12, 0, 0)))
  Sr2 <- Polygon(cbind(c(1, 1, 9, 9, 1), c(11, 12, 12, 11, 11)))
  Sr3 <- Polygon(cbind(c(1, 1, 9, 9, 1), c(0, 1, 1, 0, 0)))
  Sr5 <- Polygon(cbind(c(4, 4, 5, 5, 4), c(4, 8, 8, 4, 4)))
  
  Srs1 <- Polygons(list(Sr1), "s1")
  Srs2 <- Polygons(list(Sr2), "s2")
  Srs3 <- Polygons(list(Sr3), "s3")
  Srs4 <- Polygons(list(Sr4), "s4")
  Srs5 <- Polygons(list(Sr5), "s5")
  
  pols <- SpatialPolygons(list(Srs1, Srs2, Srs3, Srs4, Srs5), 1:5)
  xymat <- matrix(3, 3, nrow = 1, ncol = 2)
  latlonproj<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  sp::proj4string(pols) <- sp::CRS(latlonproj)
  
  expect_is(costrasterGen(xymat, pols, projstr = latlonproj), "RasterLayer")
  expect_message(costrasterGen(xymat, pols, projstr = NULL), "Warning, the projection of polygons does not match projstr. See rgdal::spTransform")
  
})

