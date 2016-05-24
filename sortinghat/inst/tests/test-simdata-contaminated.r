library('testthat')
library('sortinghat')
library('mvtnorm')

context("Generate multivariate contaminated normal populations")

test_that("Two multivariate contaminated normal populations are generated correctly", {
  seed <- 42
  
  # Generates 10 observations from each of two multivariate contaminated normal
  # populations with equal covariance matrices. Each population has a
  # contamination probability of 0.05 and scale contamination of 10.
  sample_sizes <- c(10, 10)
  means <- list(c(3, 0), c(0, 3))
  cov_identity <- diag(2)
  epsilon <- 0.05
  kappa <- 10
  
  # Data from sortinghat
  data <- simdata_contaminated(n = sample_sizes, mean = means, cov = cov_identity,
                               epsilon = epsilon, kappa = kappa, seed = seed)

  # Manually generated test data
  set.seed(seed)

  contam1 <- rbinom(sample_sizes[1], size = 1, prob = epsilon)
  normal_uncontam1 <- rmvnorm(sample_sizes[1], means[[1]], cov_identity)
  normal_contam1 <- rmvnorm(sample_sizes[1], means[[1]], kappa * cov_identity)
  x1 <- (1 - contam1) * normal_uncontam1 + contam1 * normal_contam1

  contam2 <- rbinom(sample_sizes[2], size = 1, prob = epsilon)
  normal_uncontam2 <- rmvnorm(sample_sizes[2], means[[2]], cov_identity)
  normal_contam2 <- rmvnorm(sample_sizes[2], means[[2]], kappa * cov_identity)
  x2 <- (1 - contam2) * normal_uncontam2 + contam2 * normal_contam2

  x <- rbind(x1, x2)

  y1 <- rep.int(1, times = sample_sizes[1])
  y2 <- rep.int(2, times = sample_sizes[2])
  y <- factor(c(y1, y2))
  
  # Tests that both the features and labels are equal
  expect_equal(data$x, x)
  expect_equal(data$y, y)
})

test_that("Three multivariate contaminated normal populations are generated correctly", {
  seed <- 42
  
  # Generates 10 observations from each of three multivariate contaminated
  # normal populations with unequal covariance matrices. The contamination
  # probabilities and scales differ for each population as well.
  sample_sizes <- c(10, 20, 30)
  means <- list(c = c(-3, -3), c(0, 0), c(3, 3))
  cov_identity <- diag(2)
  cov_list <- list(cov_identity, 2 * cov_identity, 3 * cov_identity)
  epsilon <- c(0.05, 0.1, 0.2)
  kappa <- c(2, 5, 10)
  
  # Data from sortinghat
  data <- simdata_contaminated(n = sample_sizes, mean = means, cov = cov_list,
                               epsilon = epsilon, kappa = kappa, seed = seed)

  # Manually generated test data
  set.seed(seed)

  contam1 <- rbinom(sample_sizes[1], size = 1, prob = epsilon[1])
  normal_uncontam1 <- rmvnorm(sample_sizes[1], means[[1]], cov_list[[1]])
  normal_contam1 <- rmvnorm(sample_sizes[1], means[[1]], kappa[1] * cov_list[[1]])
  x1 <- (1 - contam1) * normal_uncontam1 + contam1 * normal_contam1

  contam2 <- rbinom(sample_sizes[2], size = 1, prob = epsilon[2])
  normal_uncontam2 <- rmvnorm(sample_sizes[2], means[[2]], cov_list[[2]])
  normal_contam2 <- rmvnorm(sample_sizes[2], means[[2]], kappa[2] * cov_list[[2]])
  x2 <- (1 - contam2) * normal_uncontam2 + contam2 * normal_contam2

  contam3 <- rbinom(sample_sizes[3], size = 1, prob = epsilon[3])
  normal_uncontam3 <- rmvnorm(sample_sizes[3], means[[3]], cov_list[[3]])
  normal_contam3 <- rmvnorm(sample_sizes[3], means[[3]], kappa[3] * cov_list[[3]])
  x3 <- (1 - contam3) * normal_uncontam3 + contam3 * normal_contam3

  x <- rbind(x1, x2, x3)

  y1 <- rep.int(1, times = sample_sizes[1])
  y2 <- rep.int(2, times = sample_sizes[2])
  y3 <- rep.int(3, times = sample_sizes[3])
  y <- factor(c(y1, y2, y3))
  
  # Tests that both the features and labels are equal
  expect_equal(data$x, x)
  expect_equal(data$y, y)
})






