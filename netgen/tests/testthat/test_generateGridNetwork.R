context("generate grid network")

test_that("generate grid network works as expected", {
  n.points = 100L
  n.dims = 2:4
  lower = 10
  upper = 50

  checkNetwork = function(x, n.points, n.dim, lower = 0, upper = 100) {
    expect_is(x, "Network")
    expect_equal(n.points, getNumberOfNodes(x))
    expect_equal(n.dim, ncol(x$coordinates), info = paste("Number of columns does not match for n.dim", n.dim))
    expect_true(all((x$coordinates >= lower) & (x$coordinates <= upper)))
  }

  n.points.per.dim = 3L
  for (n.dim in n.dims) {
    x = generateGridNetwork(n.points.per.dim = n.points.per.dim, n.dim = n.dim, lower = lower, upper = upper)
    checkNetwork(x, n.points.per.dim^n.dim, n.dim, lower, upper)
  }
})
