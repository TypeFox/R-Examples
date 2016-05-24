context("morphing")

test_that("morphing generates a network", {
  set.seed(1)
  reps = 3L
  n.points = 30L
  n.depots = 2L
  for (i in seq(reps)) {
    x = generateRandomNetwork(n.points = n.points)
    y = generateClusteredNetwork(n.points = n.points, n.cluster = sample(2:3, 1))

    z = morphInstances(x, y, alpha = 0.5)
    expect_is(z, "Network")
    expect_equal(getNumberOfNodes(z), n.points)
    expect_true(all(z$coordinates >= 0) && all(z$coordinates <= 100))
  }

  # check that fails if we have different number of nodes
  x = generateRandomNetwork(n.points = n.points)
  y = generateRandomNetwork(n.points = 2 * n.points)
  expect_error(morphInstances(x, y, alpha = 0.5))

  # check that depots are preserved
  x = generateRandomNetwork(n.points = n.points, n.depots = n.depots)
  y = generateClusteredNetwork(n.points = n.points, n.cluster = 2L, n.depots = n.depots)
  z = morphInstances(x, y, alpha = 0.5)

  expect_is(z, "Network")
  expect_true(hasDepots(z))
  expect_equal(ncol(getDepotCoordinates(z)), n.depots)
  expect_equal(nrow(getDepotCoordinates(z)), n.depots)
  expect_equal(getNumberOfDepots(z), n.depots)

  # check user-defined point matching ...
  # ... does not work with depots
  x = generateRandomNetwork(n.points = 3L, n.depots = n.depots)
  y = generateClusteredNetwork(n.points = 3L, n.cluster = 2L, n.depots = n.depots)
  expect_error(morphInstances(x, y, alpha = 0.5, point.matching = matrix(c(1:3, 1:3), ncol = 2)))

  # ... does work with different point matching algorithms
  x = generateRandomNetwork(n.points = 5L)
  y = generateRandomNetwork(n.points = 5L)

  matching.lp = getOptimalPointMatching(x$coordinates, y$coordinates, method = "lp")
  matching.pr = getOptimalPointMatching(x$coordinates, y$coordinates, method = "push_relabel")
  expect_true(all(matching.lp == matching.pr))
  z = morphInstances(x, y, alpha = 0.5, point.matching = matching.lp)
  expect_is(z, "Network")
  z = morphInstances(x, y, alpha = 0.5, point.matching = matching.pr)
  expect_is(z, "Network")

  # ... does (not) work with faulty/correct point matching
  x = generateRandomNetwork(n.points = 5L)
  y = generateClusteredNetwork(n.points = 5L, n.cluster = 2L)
  expect_error(morphInstances(x, y, alpha = 0.5, point.matching = matrix(c(1:5, 2:6), ncol = 2)))
  z = morphInstances(x, y, alpha = 0.3, point.matching = matrix(c(1:5, sample(1:5)), ncol = 2))
  expect_is(z, "Network")

  # check point matching when passing networks
  x = generateRandomNetwork(n.points = 5L)
  matching.matrix = getOptimalPointMatching(x$coordinates, y$coordinates, method = "lp")
  matching.network = getOptimalPointMatching(x, y, method = "lp")
  expect_true(all(matching.lp == matching.pr))
})
