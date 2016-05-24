context("import and export to data formats")

# Make assertions on two networks for testing if both are equal.
# @param x, y [Network]
expect_equal_networks = function(x, y) {
  expect_equal(getNumberOfNodes(x), getNumberOfNodes(y))
  expect_equal(getNumberOfClusters(x), getNumberOfClusters(y))
  expect_true(all((x$coordinates - y$coordinates) < 0.01))
  expect_equal(x$name, y$name)
  expect_equal(x$comment, y$comment)
  if (!is.null(x$arrival.times)) {
    expect_true(all((x$arrival.times - y$arrival.times) < 0.01))
  }
}

test_that("import and export to TSPlib format is running fine", {
  # @param x Network
  testExportAndImport = function(x, name, comment) {
    x$name = name
    x$comment = comment
    fn = tempfile(fileext = ".tsp")
    exportToTSPlibFormat(x, fn, use.extended.format = FALSE)
    unlink(fn)
    exportToTSPlibFormat(x, fn, use.extended.format = TRUE)
    y = importFromTSPlibFormat(fn)
    expect_equal_networks(x, y)
  }

  x = generateRandomNetwork(n.points = 10L)
  testExportAndImport(x, name = "test", comment = "n.points=10")
  x = generateClusteredNetwork(n.points = 10L, n.cluster = 2L)
  testExportAndImport(x, name = "test2", comment = "n.points=10;n.cluster=2")

  # no name, thus we should throw an error
  x$name = NULL
  expect_error(exportToTSPlibFormat(x, "test.tsp"))

  # depots not exported at the moment
  x = generateRandomNetwork(n.points = 10L, n.depots = 2L)
  expect_error(exportToTSPlibFormat(x, "test.tsp"))
})

test_that("import and export to proprietary format is running fine", {
  # @param x Network
  testExportAndImport = function(x, name, comment) {
    x$name = name
    x$comment = comment
    fn = tempfile(fileext = ".tsp")
    exportToFile(x, fn)
    y = importFromFile(fn)
    unlink(fn)
    expect_equal_networks(x, y)
  }
  x = generateRandomNetwork(n.points = 10L)
  testExportAndImport(x, name = "test", comment = "Random network")

  x = generateClusteredNetwork(n.points = 10L, n.cluster = 2L, n.depots = 2L)
  testExportAndImport(x, name = "test", comment = "Random network with depots and clusters")

  x = dynamise(x, dyn.customers.ratio = 0.25, arrival.limit = 100)
  testExportAndImport(x, name = "test", comment = "Random network with depots and clusters and arrival times")
})
