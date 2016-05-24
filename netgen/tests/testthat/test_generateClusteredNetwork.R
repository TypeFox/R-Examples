context("generate clustered network")

test_that("generate clustered network works as expected", {
  n.cluster = 5L
  n.points = 100L
  lower = 1
  upper = 100
  n.depots = 2L

  checkClusteredInstance = function(x, n.cluster, n.points, n.depots = 0L, lower = 0, upper = 100) {
    expect_is(x, "Network")
    expect_is(x, "ClusteredNetwork")
    expect_equal(n.cluster, getNumberOfClusters(x))
    expect_true(length(setdiff(n.cluster, unique(x$membership))) == 0)
    expect_equal(n.points, getNumberOfNodes(x))
    expect_equal(n.depots, getNumberOfDepots(x))
    expect_output(print(x), "^Clustered 2-dimensional network")
    if (n.depots == 0L) {
      expect_error(getDepotCoordinates(x))
    }
    expect_true(all(x$coordinates <= upper))
    expect_true(all(x$coordinates >= lower))
  }

  # WITHOUT DEPOTS
  x = generateClusteredNetwork(n.cluster, n.points, lower = lower, upper = upper)
  checkClusteredInstance(x, n.cluster, n.points, lower = lower, upper = upper)

  # check mirroring strategy
  x = generateClusteredNetwork(n.cluster, n.points, lower = lower, upper = upper, out.of.bounds.handling = "mirror")
  checkClusteredInstance(x, n.cluster, n.points, lower = lower, upper = upper)

  #FIXME: reenable when "random.partition" is finished
  # x = generateClusteredNetwork(n.cluster, n.points, lower = lower, distribution.strategy = "random.partition")
  # checkClusteredInstance(x, n.cluster, lower = lower, n.points)

  # WITH DEPOTS
  x = generateClusteredNetwork(n.cluster, n.points, lower = lower, upper = upper, n.depots = n.depots)
  # in this case we have to nodes (the depots) more!
  checkClusteredInstance(x, n.cluster, n.points, n.depots = 2L, lower = lower, upper = upper)

  # WITH CUSTOM CLUSTER CENTERS
  x = generateClusteredNetwork(n.cluster = 10L, n.points, cluster.centers = matrix(c(10, 10, 40, 40), ncol = 2, byrow = TRUE))
  checkClusteredInstance(x, n.cluster = 2L, n.points)

  # WITH CUSTOM COVARIANCE MATRIZES
  sigma1 = 5 * diag(2)
  sigma2 = matrix(c(5, 1, 1, 5), ncol = 2, byrow = TRUE)
  sigmas = list(sigma1, sigma2)
  x = generateClusteredNetwork(n.cluster = 2L, n.points = 40L, sigmas = sigmas)
  # check plotting
  library(ggplot2)
  pl = autoplot(x)
  expect_is(pl, c("gg", "ggplot"))
})
