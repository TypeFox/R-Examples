library(survival)
library(bnnSurvival)
library(pec)

context("pec")

## Use only 1 core
options(mc.cores = 1)

## Data
dat <- veteran
n <- nrow(dat)
idx <- sample(n, 2/3*n)
train_data <- dat[idx, ]
test_data <- dat[-idx, ]
timepoints <- sort(unique(train_data$time))
formula <- formula(Surv(time, status) ~ age + trt + prior)

test_that("pec predictions works", {
  bnn <- bnnSurvival(formula, train_data, k = 5, num_base_learners = 3)
  fitpec <- pec(bnn, formula = formula, data = test_data)
  expect_that(fitpec, is_a("pec"))
})