context("dotproduct")

test_that("dotproduct functions properly with two iterators", {
  it <- iterators::iter(1:5)
  it2 <- iterators::iter(1:5)
  expect_equal(dotproduct(it, it2), 55)
})

test_that("dotproduct functions properly with two numeric vectors", {
  vec1 <- c(10, 10)
  vec2 <- c(20, 20)
  expect_equal(dotproduct(vec1, vec2), 400)
})
