context("Basic simplex in 3D")

test_that("everything is sane", {
  ordinalConstraints <- function(dim) {
    mergeConstraints(lapply(1:(dim-1), function(i) {ordinalConstraint(dim, i, i + 1)}))
  }

  n <- 3
  transform <- simplex.createTransform(n)
  constr <- simplex.createConstraints(transform, ordinalConstraints(n))
  seedPoint <- createSeedPoint(constr, homogeneous=TRUE)

  N <- 10000
  samples <- har(seedPoint, constr, N, 1, homogeneous=TRUE, transform=transform)$samples

  # Check dimension
  expect_equal(dim(samples), c(N, n))

  # Check that w_i >= w_i+1
  expect_true(all(samples[,1] >= samples[,2]))
  expect_true(all(samples[,2] >= samples[,3]))

  # Check that w_i >= 0
  expect_true(all(samples >= 0))

  # Check that sum_i w_i = 1
  expect_equal(apply(samples, 1, sum), rep(1, N))

  # Check that the points are not all identical
  expect_true(all(apply(samples, 2, sd) > 0))

  # Check that seed point is not included in sample
  expect_true(any(samples[1,] != transform %*% seedPoint))

  samples <- har(seedPoint, constr, N, 1, homogeneous=TRUE)$samples

  # Check dimension
  expect_equal(dim(samples), c(N, n))

  # Check homogeneous coordinate
  expect_equal(samples[,3], rep(1, N))

  # Check that seed point is not included in sample
  expect_true(any(samples[1,] != seedPoint))
})
