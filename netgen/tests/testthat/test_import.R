context("import of TSPlib files")

test_that("import works well for all EDGE_WEIGHT_TYPES/EDGE_WEIGHT_FORMATS", {
  testdata = system.file("testdata", package = "netgen")
  test.files = list.files(testdata, full.names = TRUE)
  for (test.file in test.files) {
    for (round.distances in c(TRUE, FALSE)) {
      x = importFromTSPlibFormat(test.file, round.distances = round.distances)
      expect_is(x, "Network", info = sprintf("Object is not a 'Network' for test file %s", basename(test.file)))
      expect_true(is.matrix(x$coordinates))
      expect_true(is.matrix(x$distance.matrix))
      expect_true(x$edge.weight.type %in% getValidEdgeWeightTypes())
    }
  }
})
