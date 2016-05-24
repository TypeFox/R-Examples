context("rescaling")

test_that("rescaling actually rescales network", {
  expect_in_bounds = function(x) {
    expect_true(all(x$coordinates >= 0) && all(x$coordinates <= 1))
    expect_true(any(x$coordinates == 0) && any(x$coordinates == 1))
    if (hasDepots(x)) {
      expect_true(all(x$depot.coordinates >= 0) && all(x$depot.coordinates <= 1))
    }
  }
  x = generateRandomNetwork(n.points = 10L)
  x.rg = range(x$coordinates)
  # GLOBAL RESCALING
  xr = rescaleNetwork(x, method = "global")

  # check for conditions that must be valid
  expect_in_bounds(xr)
  expect_true(any(xr$coordinates == 0) && any(xr$coordinates == 1))

  # scale back
  xr$coordinates = xr$coordinates * (x.rg[2] - x.rg[1]) + x.rg[1]
  expect_true(all((xr$coordinates - x$coordinates) < 0.00001))

  # GLOBAL RESCALING (2nd)
  xr = rescaleNetwork(x, method = "global2")
  expect_in_bounds(xr)

  # check for conditions that must be valid
  expect_true(any(xr$coordinates == 0) && any(xr$coordinates == 1))

  xr = rescaleNetwork(x, method = "by.dimension")
  expect_in_bounds(xr)
  expect_true(any(xr$coordinates[, 1] == 0))
  expect_true(any(xr$coordinates[, 1] == 1))
  expect_true(any(xr$coordinates[, 2] == 0))
  expect_true(any(xr$coordinates[, 2] == 1))

  # check rescaling of instances with depots
  x = generateClusteredNetwork(n.points = 50L, n.cluster = 3L, n.depots = 2L)
  xr = rescaleNetwork(x)
  expect_in_bounds(xr)
  expect_equal(getNumberOfNodes(x), getNumberOfNodes(xr))
  expect_equal(getNumberOfDepots(x), getNumberOfDepots(xr))
  expect_equal(getNumberOfClusters(x), getNumberOfClusters(xr))
})
