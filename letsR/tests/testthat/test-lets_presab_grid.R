context("Test for lets.presab.grid")

# Grid
sp.r <- rasterToPolygons(raster(resol = 5))
slot(sp.r, "data") <- cbind("ID" = 1:length(sp.r),
                            slot(sp.r, "data"))

# Species polygons
data(Phyllomedusa)
projection(Phyllomedusa) <- projection(sp.r)


test_that("lets.presab.grid works fine", {
  resu <- lets.presab.grid(Phyllomedusa, sp.r, "ID")
  expect_equal(class(resu), "list")
  expect_equal(class(resu[[1]]), "matrix")
})

