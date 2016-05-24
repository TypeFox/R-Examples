context("Basic function")

test_that("graticule creation is successful", {
  expect_that(graticule(seq(100, 240, by = 15), seq(-85, -30, by = 15)), is_a("SpatialLinesDataFrame"))
  expect_that(graticule(seq(100, 240, by = 15), seq(-85, -30, by = 15), nverts = 20), is_a("SpatialLinesDataFrame"))
  expect_that(graticule(seq(100, 240, by = 15), seq(-85, -30, by = 15), nverts = 100, xlim = c(-180, 180), ylim = c(-60, -30)), is_a("SpatialLinesDataFrame"))

})

test_that("labels work", {
  expect_that(graticule_labels(seq(100, 240, by = 15), seq(-85, -30, by = 15)), is_a("SpatialPointsDataFrame"))
})

