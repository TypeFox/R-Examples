library('testthat')
library('sortinghat')
library('mvtnorm')

context("Generate normal populations with Friedman et al.'s (1989) covariance structure")

test_that("Experiment #1 works correctly.", {
  seed <- 42
  
  # Generates observations from three multivariate normal populations having a
  # covariance structure given in Experiment #1.
  sample_sizes <- c(10, 20, 30)
  p <- 20

  # Data from sortinghat
  data <- simdata_friedman(n = sample_sizes, p = p, experiment = 1, seed = seed)

  # Manually generated test data
  cov_identity <- diag(p)
  mean1 <- rep(0, p)
  mean2 <- c(3, rep(0, p - 1))
  mean3 <- c(0, 3, rep(0, p - 2))

  set.seed(seed)
  x1 <- rmvnorm(n = sample_sizes[1], mean = mean1, sigma = cov_identity)
  x2 <- rmvnorm(n = sample_sizes[2], mean = mean2, sigma = cov_identity)
  x3 <- rmvnorm(n = sample_sizes[3], mean = mean3, sigma = cov_identity)

  x <- rbind(x1, x2, x3)

  y1 <- rep.int(1, times = sample_sizes[1])
  y2 <- rep.int(2, times = sample_sizes[2])
  y3 <- rep.int(3, times = sample_sizes[3])
  y <- factor(c(y1, y2, y3))
  
  # Tests that both the features and labels are equal
  expect_equal(data$x, x)
  expect_equal(data$y, y)
})

test_that("Experiment #2 works correctly.", {
  seed <- 42
  
  # Generates observations from three multivariate normal populations having a
  # covariance structure given in Experiment #1.
  sample_sizes <- c(15, 15, 15)
  p <- 30

  # Data from sortinghat
  data <- simdata_friedman(n = sample_sizes, p = p, experiment = 2, seed = seed)

  # Manually generated test data
  cov_identity <- diag(p)
  mean1 <- rep(0, p)
  mean2 <- c(3, rep(0, p - 1))
  mean3 <- c(0, 4, rep(0, p - 2))

  set.seed(seed)
  x1 <- rmvnorm(n = sample_sizes[1], mean = mean1, sigma = 1 * cov_identity)
  x2 <- rmvnorm(n = sample_sizes[2], mean = mean2, sigma = 2 * cov_identity)
  x3 <- rmvnorm(n = sample_sizes[3], mean = mean3, sigma = 3 * cov_identity)

  x <- rbind(x1, x2, x3)

  y1 <- rep.int(1, times = sample_sizes[1])
  y2 <- rep.int(2, times = sample_sizes[2])
  y3 <- rep.int(3, times = sample_sizes[3])
  y <- factor(c(y1, y2, y3))
  
  # Tests that both the features and labels are equal
  expect_equal(data$x, x)
  expect_equal(data$y, y)
})

test_that("Experiment #3 works correctly.", {
  seed <- 42
  
  # Generates observations from three multivariate normal populations having a
  # covariance structure given in Experiment #1.
  sample_sizes <- c(10, 20, 30)
  p <- 30

  # Data from sortinghat
  data <- simdata_friedman(n = sample_sizes, p = p, experiment = 3, seed = seed)

  # Manually generated test data
  eigenvals <- (9 * (seq_len(p) - 1) / (p - 1) + 1)^2
  cov1 <- diag(eigenvals)
  cov2 <- diag(eigenvals)
  cov3 <- diag(eigenvals)

  mean1 <- rep(0, p)
  mean2 <- 2.5 * sqrt(eigenvals / p) * (p - seq_len(p)) / (p/2 - 1)
  mean3 <- (rep(-1, p)^seq_len(p)) * mean2

  set.seed(seed)
  x1 <- rmvnorm(n = sample_sizes[1], mean = mean1, sigma = cov1)
  x2 <- rmvnorm(n = sample_sizes[2], mean = mean2, sigma = cov2)
  x3 <- rmvnorm(n = sample_sizes[3], mean = mean3, sigma = cov3)

  x <- rbind(x1, x2, x3)

  y1 <- rep.int(1, times = sample_sizes[1])
  y2 <- rep.int(2, times = sample_sizes[2])
  y3 <- rep.int(3, times = sample_sizes[3])
  y <- factor(c(y1, y2, y3))
  
  # Tests that both the features and labels are equal
  expect_equal(data$x, x)
  expect_equal(data$y, y)

  # Tests that the ratio of the max to min eigenvalues is 100
  expect_equal(max(eigenvals) / min(eigenvals), 100)
})

test_that("Experiment #4 works correctly.", {
  seed <- 42
  
  # Generates observations from three multivariate normal populations having a
  # covariance structure given in Experiment #1.
  sample_sizes <- c(10, 20, 30)
  p <- 30

  # Data from sortinghat
  data <- simdata_friedman(n = sample_sizes, p = p, experiment = 4, seed = seed)

  # Manually generated test data
  eigenvals <- (9 * (seq_len(p) - 1) / (p - 1) + 1)^2
  cov1 <- diag(eigenvals)
  cov2 <- diag(eigenvals)
  cov3 <- diag(eigenvals)

  mean1 <- rep(0, p)
  mean2 <- 2.5 * sqrt(eigenvals / p) * (seq_len(p) - 1) / (p/2 - 1)
  mean3 <- (rep(-1, p)^seq_len(p)) * mean2

  set.seed(seed)
  x1 <- rmvnorm(n = sample_sizes[1], mean = mean1, sigma = cov1)
  x2 <- rmvnorm(n = sample_sizes[2], mean = mean2, sigma = cov2)
  x3 <- rmvnorm(n = sample_sizes[3], mean = mean3, sigma = cov3)

  x <- rbind(x1, x2, x3)

  y1 <- rep.int(1, times = sample_sizes[1])
  y2 <- rep.int(2, times = sample_sizes[2])
  y3 <- rep.int(3, times = sample_sizes[3])
  y <- factor(c(y1, y2, y3))
  
  # Tests that both the features and labels are equal
  expect_equal(data$x, x)
  expect_equal(data$y, y)

  # Tests that the ratio of the max to min eigenvalues is 100
  expect_equal(max(eigenvals) / min(eigenvals), 100)
})

test_that("Experiment #5 works correctly.", {
  seed <- 42
  
  # Generates observations from three multivariate normal populations having a
  # covariance structure given in Experiment #1.
  sample_sizes <- c(10, 20, 30)
  p <- 30

  # Data from sortinghat
  data <- simdata_friedman(n = sample_sizes, p = p, experiment = 5, seed = seed)

  # Manually generated test data
  eigenvals1 <- (9 * (seq_len(p) - 1) / (p - 1) + 1)^2
  eigenvals2 <- (9 * (p - seq_len(p)) / (p - 1) + 1)^2
  eigenvals3 <- (9 * (seq_len(p) - (p - 1) / 2) / (p - 1))^2
  
  cov1 <- diag(eigenvals1)
  cov2 <- diag(eigenvals2)
  cov3 <- diag(eigenvals3)

  mean1 <- rep(0, p)
  mean2 <- rep(0, p)
  mean3 <- rep(0, p)

  set.seed(seed)
  x1 <- rmvnorm(n = sample_sizes[1], mean = mean1, sigma = cov1)
  x2 <- rmvnorm(n = sample_sizes[2], mean = mean2, sigma = cov2)
  x3 <- rmvnorm(n = sample_sizes[3], mean = mean3, sigma = cov3)

  x <- rbind(x1, x2, x3)

  y1 <- rep.int(1, times = sample_sizes[1])
  y2 <- rep.int(2, times = sample_sizes[2])
  y3 <- rep.int(3, times = sample_sizes[3])
  y <- factor(c(y1, y2, y3))
  
  # Tests that both the features and labels are equal
  expect_equal(data$x, x)
  expect_equal(data$y, y)
})

test_that("Experiment #6 works correctly.", {
  seed <- 42
  
  # Generates observations from three multivariate normal populations having a
  # covariance structure given in Experiment #1.
  sample_sizes <- c(10, 20, 30)
  p <- 30

  # Data from sortinghat
  data <- simdata_friedman(n = sample_sizes, p = p, experiment = 6, seed = seed)

  # Manually generated test data
  eigenvals1 <- (9 * (seq_len(p) - 1) / (p - 1) + 1)^2
  eigenvals2 <- (9 * (p - seq_len(p)) / (p - 1) + 1)^2
  eigenvals3 <- (9 * (seq_len(p) - (p - 1) / 2) / (p - 1))^2
  
  cov1 <- diag(eigenvals1)
  cov2 <- diag(eigenvals2)
  cov3 <- diag(eigenvals3)

  mean1 <- rep(0, p)
  mean2 <- 14 / sqrt(p) * rep(1, p)
  mean3 <- (rep(-1, p)^seq_len(p)) * mean2

  set.seed(seed)
  x1 <- rmvnorm(n = sample_sizes[1], mean = mean1, sigma = cov1)
  x2 <- rmvnorm(n = sample_sizes[2], mean = mean2, sigma = cov2)
  x3 <- rmvnorm(n = sample_sizes[3], mean = mean3, sigma = cov3)

  x <- rbind(x1, x2, x3)

  y1 <- rep.int(1, times = sample_sizes[1])
  y2 <- rep.int(2, times = sample_sizes[2])
  y3 <- rep.int(3, times = sample_sizes[3])
  y <- factor(c(y1, y2, y3))
  
  # Tests that both the features and labels are equal
  expect_equal(data$x, x)
  expect_equal(data$y, y)
})

