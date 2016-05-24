library('testthat')
library('sortinghat')
library('mvtnorm')

context("Generate multivariate normal populations")

test_that("Two multivariate normal populations are generated correctly", {
  seed <- 42
  
  # Generates observations from each of two multivariate normal populations with
  # equal covariance matrices.
  sample_sizes <- c(10, 10)
  means <- list(c(1, 0), c(0, 1))
  cov_identity <- diag(2)

  # Data from sortinghat
  data <- simdata_normal(n = sample_sizes, mean = means, cov = cov_identity,
                         seed = seed)

  # Manually generated test data
  set.seed(seed)
  x1 <- rmvnorm(n = sample_sizes[1], mean = means[[1]], sigma = cov_identity)
  x2 <- rmvnorm(n = sample_sizes[2], mean = means[[2]], sigma = cov_identity)
  x <- rbind(x1, x2)

  y1 <- rep.int(1, times = sample_sizes[1])
  y2 <- rep.int(2, times = sample_sizes[2])
  y <- factor(c(y1, y2))
  
  # Tests that both the features and labels are equal
  expect_equal(data$x, x)
  expect_equal(data$y, y)
})

test_that("Three multivariate normal populations are generated correctly", {
  seed <- 42

  # Generates observations from each of three multivariate normal populations
  # with unequal covariance matrices.
  sample_sizes <- c(10, 20, 30)
  means <- list(c = c(-3, -3), c(0, 0), c(3, 3))
  cov_identity <- diag(2)
  cov_list <- list(cov_identity, 2 * cov_identity, 3 * cov_identity)

  # Data from sortinghat
  data <- simdata_normal(n = sample_sizes, mean = means, cov = cov_list,
                         seed = seed)

  # Manually generated test data
  set.seed(seed)
  x1 <- rmvnorm(n = sample_sizes[1], mean = means[[1]], sigma = cov_list[[1]])
  x2 <- rmvnorm(n = sample_sizes[2], mean = means[[2]], sigma = cov_list[[2]])
  x3 <- rmvnorm(n = sample_sizes[3], mean = means[[3]], sigma = cov_list[[3]])
  x <- rbind(x1, x2, x3)

  y1 <- rep.int(1, times = sample_sizes[1])
  y2 <- rep.int(2, times = sample_sizes[2])
  y3 <- rep.int(3, times = sample_sizes[3])
  y <- factor(c(y1, y2, y3))
  
  # Tests that both the features and labels are equal
  expect_equal(data$x, x)
  expect_equal(data$y, y)
})
