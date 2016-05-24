context("generate random instance")

test_that("generate random instance works as expected", {
  n.points = 100L
  n.dims = 2:4
  for (n.dim in n.dims) {
    x = generateRandomNetwork(n.points, n.dim = n.dim)
    expect_is(x, "Network")
    expect_equal(n.points, getNumberOfNodes(x))
    expect_equal(n.dim, ncol(x$coordinates), info = paste("Number of columns does not match for n.dim", n.dim))
  }

  expect_error(generateRandomNetwork(n.points = 10, lower = 10, upper = 10))

  # check if points are within bounds
  lower = 0
  upper = 0.5

  x = generateRandomNetwork(n.points = 10L, lower = lower, upper = upper)
  expect_true(all((x$coordinates >= lower) & (x$coordinates <= upper)))

  x = generateRandomNetwork(n.points = 10L, n.depots = 1L, lower = lower, upper = upper)
  expect_true(all((x$coordinates >= lower) & (x$coordinates <= upper)))
  expect_equal(1, getNumberOfDepots(x))
  expect_true(isEuclidean(x))
})
