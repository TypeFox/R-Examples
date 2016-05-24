library('testthat')
library('sortinghat')

context("Partitioning a Data Set for Cross-validation")

test_that("Defaults work correctly on Fisher's Iris data set.", {
  class_labels <- iris$Species
  num_folds <- 10
  n <- length(class_labels)

  hold_out <- n / num_folds

  set.seed(42)
  cv_out <- cv_partition(iris$Species)

  # Tests that the number of training and test observations in each fold
  # is correct. This assumes that the number of folds is 10 by default.
  num_training <- sapply(lapply(cv_out, "[[", "training"), length)
  num_test <- sapply(lapply(cv_out, "[[", "test"), length)
  
  expect_true(all(num_training == n - hold_out))
  expect_true(all(num_test == hold_out))
})

test_that("Specifying num_folds works correctly on Fisher's Iris data set.", {
  class_labels <- iris$Species
  num_folds <- 5
  n <- length(class_labels)

  hold_out <- n / num_folds
  
  cv_out <- cv_partition(y = class_labels, num_folds = 5, seed = 42)

  # Tests that the number of training and test observations in each fold
  # is correct.
  num_training <- sapply(lapply(cv_out, "[[", "training"), length)
  num_test <- sapply(lapply(cv_out, "[[", "test"), length)
  
  expect_true(all(num_training == n - hold_out))
  expect_true(all(num_test == hold_out))
})

test_that("Specifying hold_out works correctly on Fisher's Iris data set.", {
  class_labels <- iris$Species
  hold_out <- 15
  n <- length(class_labels)
  
  cv_out <- cv_partition(y = class_labels, hold_out = hold_out, seed = 42)

  # Tests that the number of training and test observations in each fold
  # is correct.
  num_training <- sapply(lapply(cv_out, "[[", "training"), length)
  num_test <- sapply(lapply(cv_out, "[[", "test"), length)
  
  expect_true(all(num_training == n - hold_out))
  expect_true(all(num_test == hold_out))
})

test_that("When hold_out is provided, then num_folds is ignored.", {
  class_labels <- iris$Species
  hold_out <- 15
  n <- length(class_labels)
  
  cv_out <- cv_partition(y = class_labels, num_folds = 5, hold_out = hold_out, seed = 42)

  # Tests that the number of training and test observations in each fold
  # is correct.
  num_training <- sapply(lapply(cv_out, "[[", "training"), length)
  num_test <- sapply(lapply(cv_out, "[[", "test"), length)
  
  expect_true(all(num_training == n - hold_out))
  expect_true(all(num_test == hold_out))
})

test_that("Folds are correct when length(y) and num_folds are relatively prime", {
  class_labels <- iris$Species
  num_folds <- 17
  n <- length(class_labels)

  hold_out <- ceiling(n / num_folds)

  remainder <- n %% num_folds

  # The length of each fold should like these
  folds_test <- c(rep(hold_out, remainder), rep(hold_out - 1, num_folds - remainder))
  folds_training <- n - folds_test

  cv_out <- cv_partition(iris$Species, num_folds = num_folds, seed = 42)

  # Tests that the number of training and test observations in each fold
  # is correct. This assumes that the number of folds is 10 by default.
  num_training <- sapply(lapply(cv_out, "[[", "training"), length)
  num_test <- sapply(lapply(cv_out, "[[", "test"), length)
  
  expect_true(all(num_training == folds_training))
  expect_true(all(num_test == folds_test))
})

# Problem Case - Issue #8
test_that("No empty test folds are given when num_folds exceeds length(y)", {
  y <- structure(c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
                 .Label = c("2", "12"),
                 class = "factor")
  cv_out <- cv_partition(y, num_folds = 10, seed = 42)

  # Tests that the number of training and test observations in each fold
  # is positive
  num_training <- sapply(lapply(cv_out, "[[", "training"), length)
  num_test <- sapply(lapply(cv_out, "[[", "test"), length)
  
  expect_true(all(num_training > 0))
  expect_true(all(num_test > 0))
})
