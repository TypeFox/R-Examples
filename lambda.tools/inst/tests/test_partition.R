# :vim set filetype=R

context("partition")
test_that("Basic functionality of partition works for a small sequence", {
  x <- 1:5
  y <- matrix(c(1, 3, 5, 7, 5, 7, 9, 5), ncol=2)
  anynames(y) <- c("left", "right")
  expect_equal(partition(x, metric=sum, radius=2), y)
})

test_that("Partition works for a non-default radius", {
  x <- 1:5
  y <- matrix(c(1, 3, 6, 9, 9, 12, 9, 5), ncol=2)
  anynames(y) <- c("left", "right")
  expect_equal(partition(x, metric=sum, radius=3), y)
})

test_that("A matrix is not allowed as input", {
  x <- matrix(rnorm(10), ncol=2)
  expect_error(partition(x), "No valid function for")
})

