library('testthat')
library('sortinghat')
library('mvtnorm')

context("Generate multivariate t populations")

test_that("Two multivariate t populations are generated correctly", {
  seed <- 42
  
  # Generates observations from each of two multivariate t populations with
  # equal covariance matrices and equal degrees of freedom.
  sample_sizes <- c(10, 10)
  centroids <- list(c(3, 0), c(0, 3))
  cov_identity <- diag(2)
  df <- 4
  
  # Data from sortinghat
  data <- simdata_t(n = sample_sizes, centroid = centroids, cov = cov_identity,
                    df = df, seed = seed)

  # Manually generated test data
  set.seed(seed)
  x1 <- rmvt(n = sample_sizes[1], delta = centroids[[1]], sigma = cov_identity,
             df = df)
  x2 <- rmvt(n = sample_sizes[2], delta = centroids[[2]], sigma = cov_identity,
             df = df)
  x <- rbind(x1, x2)

  y1 <- rep.int(1, times = sample_sizes[1])
  y2 <- rep.int(2, times = sample_sizes[2])
  y <- factor(c(y1, y2))
  
  # Tests that both the features and labels are equal
  expect_equal(data$x, x)
  expect_equal(data$y, y)
})

test_that("Three multivariate t populations are generated correctly", {
  seed <- 42
  
  # Generates observations from each of three multivariate t populations
  # with unequal covariance matrices and unequal degrees of freedom.
  sample_sizes <- c(10, 20, 30)
  centroids <- list(c = c(-3, -3), c(0, 0), c(3, 3))
  cov_identity <- diag(2)
  cov_list <- list(cov_identity, 2 * cov_identity, 3 * cov_identity)
  df <- c(3, 4, 5)

  # Data from sortinghat
  data <- simdata_t(n = sample_sizes, centroid = centroids, cov = cov_list,
                    df = df, seed = seed)

  # Manually generated test data
  set.seed(seed)
  x1 <- rmvt(n = sample_sizes[1], delta = centroids[[1]], sigma = cov_list[[1]],
             df = df[1])
  x2 <- rmvt(n = sample_sizes[2], delta = centroids[[2]], sigma = cov_list[[2]],
             df = df[2])
  x3 <- rmvt(n = sample_sizes[3], delta = centroids[[3]], sigma = cov_list[[3]],
             df = df[3])
  x <- rbind(x1, x2, x3)

  y1 <- rep.int(1, times = sample_sizes[1])
  y2 <- rep.int(2, times = sample_sizes[2])
  y3 <- rep.int(3, times = sample_sizes[3])
  y <- factor(c(y1, y2, y3))
  
  # Tests that both the features and labels are equal
  expect_equal(data$x, x)
  expect_equal(data$y, y)
})






