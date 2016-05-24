context("optimizeSubInts")

test_that("optimizeSubInts", {

  f = function(x) sin(x) * x
  z = optimizeSubInts(f, interval = c(0, 50), nsub = 200L)
  fopt = f(pi * 3 / 2 + 14 * pi)
  expect_true(abs(fopt -  z$objective) < 1e-1)

  # test with nsub = 1, had a bug here
  f = function(x) sum(x^2)
  z = optimizeSubInts(f, interval = c(-10, 10), nsub = 1)
  expect_true(abs(z$minimum) < 1e-5)
})

