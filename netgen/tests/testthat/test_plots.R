context("plots")

test_that("all of our plots produce nice gglot2 objects", {
  # helper function
  expect_is_ggplot = function(pl) {
    expect_is(pl, c("gg", "ggplot"))
  }

  # more than three dimensions not "plotable"
  x = generateRandomNetwork(n.points = 10L, n.dim = 3L)
  expect_error(autoplot(x))

  # no cluster, no depots
  x = generateRandomNetwork(n.points = 20L)
  expect_is_ggplot(autoplot(x))
  expect_is_ggplot(autoplot(x, path = sample(1:20)))
  expect_is_ggplot(autoplot(x, path = sample(1:10), close.path = TRUE, path.colour = "tomato"))

  # no cluster, with depots
  x = generateRandomNetwork(n.points = 20L, n.depots = 2L, upper = 100)
  expect_is_ggplot(autoplot(x))

  # some clusters, no depots
  x = generateClusteredNetwork(n.points = 20L, n.cluster = 2L)
  expect_is_ggplot(autoplot(x))

  # some clusters, with depots
  x = generateClusteredNetwork(n.points = 20L, n.cluster = 2L, n.depots = 2L)
  expect_is_ggplot(autoplot(x))

  # autoplot morphed instance
  x = generateRandomNetwork(n.points = 10L, upper = 20)
  y = generateClusteredNetwork(n.points = 10L, n.cluster = 2L, upper = 60)
  z = morphInstances(x, y, alpha = 0.3)
  expect_is_ggplot(autoplot(z))

  for (show.arrows in c(TRUE, FALSE)) {
    for (do.facets in c(TRUE, FALSE)) {
      pl = visualizeMorphing(x, y, arrows = show.arrows, in.one.plot = do.facets)
      expect_is_ggplot(pl)
    }
  }

  # test visualization of point matchings
  x = generateRandomNetwork(n.points = 10L)
  y = generateClusteredNetwork(n.points = 10L, n.cluster = 2L)
  pm = netgen::getOptimalPointMatching(x$coordinates, y$coordinates)
  pl = visualizePointMatching(x, y, pm, highlight.longest = 3L)
  expect_is_ggplot(pl)
})
