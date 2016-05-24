context("tsearchn")
test_that("tsearchn gives the expected output", {
  x <- c(-1, -1, 1)
  y <- c(-1, 1, -1)
  p <- cbind(x, y)
  tri <- matrix(c(1, 2, 3), 1, 3)
  ## Should be in triangle #1
  ts <- tsearchn(p, tri, cbind(-1, -1),fast=FALSE)
  expect_that(ts$idx, equals(1))
  expect_that(ts$p, equals(cbind(1, 0, 0)))
  ## Should be in triangle #1
  ts <- tsearchn(p, tri, cbind(1, -1), fast=FALSE)
  expect_that(ts$idx, equals(1))
  expect_that(ts$p, equals(cbind(0, 0, 1)))
  ## Should be in triangle #1
  ts <- tsearchn(p, tri, cbind(-1, 1), fast=FALSE)
  expect_that(ts$idx, equals(1))
  expect_that(ts$p, equals(cbind(0, 1, 0)))
  ## Centroid
  ts <- tsearchn(p, tri, cbind(-1/3, -1/3), fast=FALSE)
  expect_that(ts$idx, equals(1))
  expect_that(ts$p, equals(cbind(1/3, 1/3, 1/3)))
  ## Should be outside triangle #1, so should return NA
  ts <- tsearchn(p, tri, cbind(1, 1), fast=FALSE)
  expect_true(is.na(ts$idx))
  expect_true(all(is.na(ts$p)))

  ## Create a mesh with a zero-area element (degenerate simplex)
  p <- cbind(c(-1, -1, 0, 1, 2),
             c(-1,  1, 0, 0, 0))
  tri <- rbind(c(1, 2, 3),
               c(3, 4, 5))
  ## Look for one point in one of the simplices and a point outwith the
  ## simplices. This forces tsearchn to look in all simplices. It
  ## shouldn't fail on the degenerate simplex.
  expect_warning(ts <- tsearchn(p, tri, rbind(c(-0.5, 0), c(3, 1)), fast=FALSE))
  expect_equal(ts$idx, c(1, NA))
  ts <- tsearchn(p, tri, rbind(c(-0.5, 0), c(3, 1)), fast=TRUE)
  expect_equal(ts$idx, c(1, NA))

})



