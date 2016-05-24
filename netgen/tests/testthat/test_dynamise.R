context("dynamize instances")

test_that("dynamize works as expected", {
  # setup
  n.points = 100L
  n.dynamic = 20L
  dyn.customers.ratios = c(0.25, 0.5, 0.75)
  arrival.limit = 20

  # check n.dynamic
  x = generateRandomNetwork(n.points = 100L)
  expect_error(dynamise(x))
  x = dynamise(x, n.dynamic = n.dynamic, arrival.limit = arrival.limit)
  expect_false(is.null(x$arrival.times))
  expect_equal(sum(x$arrival.times > 0), n.dynamic)

  # check dyn.customers.ratio
  for (dyn.customers.ratio in dyn.customers.ratios) {
    x = generateRandomNetwork(n.points = n.points, n.depots = 2L)
    x = dynamise(x, dyn.customers.ratio = dyn.customers.ratio, arrival.limit = arrival.limit)
    expect_false(is.null(x$arrival.times))
    expected.n.dynamic = dyn.customers.ratio * n.points
    current.n.dynamic = sum(x$arrival.times > 0)
    expect_equal(current.n.dynamic, expected.n.dynamic,
      info = sprintf("Number of dyanmic customers wrong! IS: %i, EXPECTED: %i",
        current.n.dynamic, expected.n.dynamic))
  }
})
