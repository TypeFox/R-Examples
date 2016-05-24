context("quantify an iterator of logical values")

test_that("quantify functions properly with a logical vector", {
  set.seed(42)
  x <- sample(c(TRUE, FALSE), size=10, replace=TRUE)
  expect_equal(quantify(x), sum(x))
})

test_that("quantify functions properly with an iterator from a logical vector", {
  set.seed(42)
  x <- sample(c(TRUE, FALSE), size=10, replace=TRUE)
  it <- iterators::iter(x)
  expect_equal(quantify(it), sum(x))
})
