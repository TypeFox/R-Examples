context("Checking Karcher mean")
test_that("Karcher mean with zero-row gives correct output ", {
  expect_that(karcher.mean.sphere(matrix(NA, nrow=0, ncol=2)), equals(c(phi=NA, lambda=NA)))
  expect_that(karcher.mean.sphere(matrix(NA, nrow=0, ncol=2), var=TRUE), equals(list(mean=c(phi=NA, lambda=NA), var=c(phi=NA, lambda=NA))))
})

test_that("Karcher mean with one-row gives correct output ", {
  x <- cbind(phi=1, lambda=0)
  expect_that(karcher.mean.sphere(x), equals(c(phi=1, lambda=0)))
  expect_that(karcher.mean.sphere(x, var=TRUE), equals(list(mean=c(phi=1, lambda=0), var=c(phi=NA, lambda=NA))))
})

test_that("Karcher mean with two-rows gives correct output ", {
  x <- cbind(phi=c(0, 0), lambda=c(0, pi/2))
  expect_that(karcher.mean.sphere(x), equals(c(phi=0, lambda=pi/4)))
  expect_that(karcher.mean.sphere(x, var=TRUE), equals(list(mean=c(phi=0, lambda=pi/4), var=(pi/4)^2)))
})
