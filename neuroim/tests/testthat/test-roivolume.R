
test_that("can create a RegionSphere", {
  sp1 <- BrainSpace(c(10,10,10), c(1,1,1))
  sphere <- RegionSphere(sp1, c(5,5,5), 3)
  expect_that(sphere, is_a("ROIVolume"))
})

test_that("can create a RegionCube", {
  sp1 <- BrainSpace(c(10,10,10), c(1,1,1))
  cube <- RegionCube(sp1, c(5,5,5), 3)
  expect_that(cube, is_a("ROIVolume"))
})

test_that("can convert region coordinates to indices", {
  sp1 <- BrainSpace(c(10,10,10), c(1,1,1))
  cube <- RegionCube(sp1, c(5,5,5), 3)
  idx <- gridToIndex(space(cube), coords(cube))
  vox <- indexToGrid(space(cube), idx)
  expect_equivalent(vox, coords(cube))
})

