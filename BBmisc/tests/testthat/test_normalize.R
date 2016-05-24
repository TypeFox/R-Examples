context("normalize")

test_that("normalize", {
  # vector
  x = runif(20)
  y = normalize(x, method = "range")
  expect_is(y, "numeric")
  expect_equal(range(y), c(0, 1))
  y = normalize(x, method = "range", range = c(-4, 2))
  expect_is(y, "numeric")
  expect_equal(range(y), c(-4, 2))
  y = normalize(x, method = "center")
  expect_is(y, "numeric")
  expect_equal(mean(y), 0)
  y = normalize(x, method = "standardize")
  expect_is(y, "numeric")
  expect_equal(mean(y), 0)
  expect_equal(sd(y), 1)

  # matrix
  x = matrix(runif(100), nrow = 5)
  y = normalize(x, margin = 1L)
  expect_is(y, "matrix")
  apply(y, 1, function(v) expect_equal(mean(v), 0))
  apply(y, 1, function(v) expect_equal(sd(v), 1))
  y = normalize(x, margin = 2L)
  apply(y, 2, function(v) expect_equal(mean(v), 0))
  apply(y, 2, function(v) expect_equal(sd(v), 1))

  # data.frame
  y = normalize(iris, method = "range", range = c(3, 4))
  expect_is(y, "data.frame")
  for (i in 1:4)
    expect_equal(range(y[, i]), c(3, 4))
  y[, 5L] = iris$Specis

  # constant vectors
  x = rep(1, 10)
  y = normalize(x, method = "center", on.constant = "quiet")
  expect_is(y, "numeric")
  expect_equal(y, x - x)
  y = normalize(x, method = "scale", on.constant = "quiet")
  expect_is(y, "numeric")
  expect_equal(y, x)
  y = normalize(x, method = "standardize", on.constant = "quiet")
  expect_is(y, "numeric")
  expect_equal(y, x - x)
  y = normalize(x, method = "range", on.constant = "quiet", range = c(-3, 2))
  expect_is(y, "numeric")
  expect_equal(y, rep(-0.5, 10))
  expect_error(normalize(x, method = "center", on.constant = "stop"))
  expect_error(normalize(x, method = "scale", on.constant = "stop"))
  expect_error(normalize(x, method = "standardize", on.constant = "stop"))
  expect_error(normalize(x, method = "range", on.constant = "stop"))
  expect_warning(normalize(x, method = "center", on.constant = "warn"))
  expect_warning(normalize(x, method = "scale", on.constant = "warn"))
  expect_warning(normalize(x, method = "standardize", on.constant = "warn"))
  expect_warning(normalize(x, method = "range", on.constant = "warn"))
})

test_that("normalize works with NAs", {
  # vector
  x = c(1, 2, NA)
  y = normalize(x, method = "range")
  expect_equal(y, c(0, 1, NA))
  y = normalize(x, method = "center")
  expect_equal(y, c(-0.5, 0.5, NA))

  # matrix
  x = matrix(c(1, 2, 1, NA), nrow = 2L)
  y = normalize(x, margin = 2L, method = "range")
  expect_equal(y, matrix(c(0, 1, 0.5, NA), nrow = 2L))
})
