context("as.data.frame, as.matrix, ...")

test_that("convertion works as expected", {
  checkDataFrame = function(x, n.expected.points, n.expected.dim = 2L, n.expected.depots = 0L, n.extra.cols = 1L) {
    df = as.data.frame(x)
    expect_equal(nrow(df), n.expected.points + n.expected.depots)
    expect_equal(ncol(df), n.expected.dim + n.extra.cols)

    df = as.data.frame(x, include.extra = FALSE)
    expect_equal(nrow(df), n.expected.points + n.expected.depots)
    expect_equal(ncol(df), n.expected.dim)

    mat = as.matrix(x)
    expect_equal(nrow(mat), n.expected.points + n.expected.depots)
    expect_equal(ncol(mat), n.expected.dim)
  }

  x = generateRandomNetwork(n.points = 10L)
  checkDataFrame(x, 10L)

  x = generateRandomNetwork(n.points = 10L, n.depots = 2L)
  checkDataFrame(x, 10L, n.expected.depots = 2L)

  x = generateClusteredNetwork(n.points = 10L, n.cluster = 2L)
  checkDataFrame(x, 10L, n.extra.cols = 2L)

  x = generateClusteredNetwork(n.points = 10L, n.cluster = 2L, n.depots = 2L)
  checkDataFrame(x, 10L, n.extra.cols = 2L, n.expected.depots = 2L)
})
