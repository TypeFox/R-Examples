context("pathdistGen")

test_that("pathdistGen works", {
	skip_on_cran()
	skip_on_travis()
	skip_on_appveyor()
  
  expect_error(pathdistGen(1, raster::raster(1), 1, progressbar = FALSE), "spdf object must be of class SpatialPointsDataFrame")
  
})