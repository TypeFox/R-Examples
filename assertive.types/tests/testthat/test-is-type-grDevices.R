test_that("test.is_raster.a_matrix.returns_false", {
  m <- matrix(hcl(0, 80, seq(50, 80, 10)), nrow = 4, ncol = 5)
  expect_false(is_raster(m))
})

test_that("test.is_raster.a_raster.returns_true", {
  m <- matrix(hcl(0, 80, seq(50, 80, 10)), nrow = 4, ncol = 5)
  expect_true(is_raster(as.raster(m)))
})
