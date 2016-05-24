# test_Ckmeans.1d.dp.R
#
# Joe Song
# Created: May 3, 2016

library(testthat)
library(Ckmeans.1d.dp)
context("Checking on several examples")


test_that("Testing Ckmeans.1d.dp()", {

  x <- c(-2.5, -2.5, -2.5, -2.5)
  res <- Ckmeans.1d.dp(x, 1)

  cluster.truth <- c(1, 1, 1, 1)
  expect_equal(res$cluster, cluster.truth)
  centers.truth <- c(-2.5)
  expect_equal(res$centers, centers.truth)
  withinss.truth <- c(0)
  expect_equal(res$withinss, withinss.truth)
  size.truth <- c(4)
  expect_equal(res$size, size.truth)

  x <- c(3, 3, 3, 3, 1, 1, 1, 2, 2, 2)
  res <- Ckmeans.1d.dp(x, 3)

  cluster.truth <- c(3, 3, 3, 3, 1, 1, 1, 2, 2, 2)
  expect_equal(res$cluster, cluster.truth)
  centers.truth <- c(1, 2, 3)
  expect_equal(res$centers, centers.truth)
  withinss.truth <- c(0, 0, 0)
  expect_equal(res$withinss, withinss.truth)
  size.truth <- c(3, 3, 4)
  expect_equal(res$size, size.truth)

  x <- c(3.5, 3.6, 3.7, 3.1, 1.1, 0.9, 0.8, 2.2, 1.9, 2.1)
  res <- Ckmeans.1d.dp(x, k=c(2,5))

  cluster.truth <- c(3, 3, 3, 3, 1, 1, 1, 2, 2, 2)
  expect_equal(res$cluster, cluster.truth)
  centers.truth <- c(0.933333333333, 2.066666666667, 3.475000000000)
  expect_equal(res$centers, centers.truth)
  withinss.truth <- c(0.0466666666667, 0.0466666666667, 0.2075000000000)
  expect_equal(res$withinss, withinss.truth)
  size.truth <- c(3, 3, 4)
  expect_equal(res$size, size.truth)

  x <- c(-3, 2.2, -6, 7, 9, 11, -6.3, 75, 82.6, 32.3, -9.5, 62.5, 7, 95.2)
  res <- Ckmeans.1d.dp(x, k=8)

  cluster.truth <- c(2, 2, 1, 3, 3, 3, 1, 6, 7, 4, 1, 5, 3, 8)
  expect_equal(res$cluster, cluster.truth)
  centers.truth <- c(-7.266666667, -0.400000000, 8.500000000, 32.300000000,
                     62.500000000, 75.000000000, 82.600000000, 95.200000000)
  expect_equal(res$centers, centers.truth)
  withinss.truth <- c(7.526666667, 13.520000000, 11.000000000, 0.000000000,
                      0.000000000, 0.000000000, 0.000000000, 0.000000000)
  expect_equal(res$withinss, withinss.truth)
  size.truth <- c(3, 2, 4, 1, 1, 1, 1, 1)
  expect_equal(res$size, size.truth)

  x <- cos((-10:10))
  res <- Ckmeans.1d.dp(x)
  # format(res, digits=10)

  cluster.truth <- c(1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1)
  expect_equal(res$cluster, cluster.truth)
  centers.truth <- c(-0.6592474631, 0.6751193405)
  expect_equal(res$centers, centers.truth)
  withinss.truth <- c(1.0564793100, 0.6232976959)
  expect_equal(res$withinss, withinss.truth)
  size.truth <- c(12, 9)
  expect_equal(res$size, size.truth)

  x <- dgamma(seq(1,10, by=0.5), shape=2, rate=1)
  res <- Ckmeans.1d.dp(x)
  # format(res, digits=10)

  cluster.truth <- c(3, 3, 3, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  expect_equal(res$cluster, cluster.truth)
  centers.truth <- c(0.01702193495, 0.15342151455, 0.32441508262)
  expect_equal(res$centers, centers.truth)
  withinss.truth <- c(0.006126754998, 0.004977009034, 0.004883305120)
  expect_equal(res$withinss, withinss.truth)
  size.truth <- c(13, 3, 3)
  expect_equal(res$size, size.truth)

})
