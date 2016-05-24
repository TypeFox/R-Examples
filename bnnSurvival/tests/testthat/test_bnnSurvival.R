library(survival)
library(bnnSurvival)

context("bnnSurvival")

## Use only 1 core
options(mc.cores = 1)

## Data
dat <- veteran
n <- nrow(dat)
idx <- sample(n, 2/3*n)
train_data <- dat[idx, ]
test_data <- dat[-idx, ]
timepoints <- sort(unique(train_data$time))
formula <- formula(Surv(time, status) ~ .)

## Run bnn Survival method
bnn <- bnnSurvival(formula, train_data, k = 5, num_base_learners = 3)
pred_bnn <- predict(bnn, test_data)

## Run knn Survival method
knn <- bnnSurvival(formula, train_data, k = 5, num_base_learners = 1, replace = FALSE, sample_fraction = 1)
pred_knn <- predict(knn, test_data)

test_that("result is of class bnnSurvival", {
  expect_that(bnn, is_a("bnnSurvivalEnsemble"))
  expect_that(knn, is_a("bnnSurvivalEnsemble"))
})

test_that("prediction is of correct size", {
  expect_that(dim(predictions(pred_bnn)), equals(c(nrow(test_data), length(timepoints))))
  expect_that(dim(predictions(pred_knn)), equals(c(nrow(test_data), length(timepoints))))
})







