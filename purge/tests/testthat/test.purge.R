library(purge)

context("Purged model predictions")

purge_test_helper <- function(unpurged.model, purged.model, test.new.data,
                              predict.method=predict) {
  preds.purged <- predict.method(purged.model, test.new.data)
  expect_equal(length(preds.purged), nrow(test.new.data))
  preds.unpurged <- predict.method(unpurged.model, test.new.data)
  expect_equal(preds.purged, preds.unpurged)
}

test_that("lm purge works correctly", {
  sample.size <- 1000
  x <- rnorm(sample.size)
  y <- rnorm(sample.size)
  unpurged.model <- lm(y ~ x)
  purged.model <- purge(unpurged.model)
  test.new.data <- data.frame(x=1:10)
  expect_is(purged.model, 'lm')
  purge_test_helper(unpurged.model, purged.model, test.new.data)
})

test_that("glm purge works correctly", {
  sample.size <- 1000
  x <- rnorm(sample.size)
  y <- as.factor(runif(sample.size) > 0.5)
  unpurged.model <- glm(y ~ x, family=binomial())
  purged.model <- purge(unpurged.model)
  test.new.data <- data.frame(x=1:10)
  expect_is(purged.model, 'glm')
  purge_test_helper(unpurged.model, purged.model, test.new.data)
})

test_that("merMod purge works correctly", {
  if (requireNamespace('lme4', quietly=TRUE)) {
    sample.size <- 1000
    x <- rnorm(sample.size)
    y <- rnorm(sample.size)
    z <- as.factor(runif(sample.size) > 0.5)
    unpurged.model <- lme4::lmer(y ~ x + (1|z))
    purged.model <- purge(unpurged.model)
    test.new.data <- data.frame(x=1:10, z=as.factor(runif(10) > 0.5))
    expect_is(purged.model, 'merMod')
    purge_test_helper(unpurged.model, purged.model, test.new.data)
  }
})

test_that("glmerMod purge works correctly", {
  if (requireNamespace('lme4', quietly=TRUE)) {
    sample.size <- 1000
    x <- rnorm(sample.size)
    y <- as.factor(runif(sample.size) > 0.5)
    z <- as.factor(runif(sample.size) > 0.5)
    unpurged.model <- lme4::glmer(y ~ x + (1|z), family=binomial())
    purged.model <- purge(unpurged.model)
    test.new.data <- data.frame(x=1:10, z=as.factor(runif(10) > 0.5))
    expect_is(purged.model, 'glmerMod')
    purge_test_helper(unpurged.model, purged.model, test.new.data)
  }
})

test_that("rpart purge works correctly", {
  if (requireNamespace('rpart', quietly=TRUE)) {
    sample.size <- 1000
    x <- rnorm(sample.size)
    y <- x + rnorm(sample.size)
    unpurged.model <- rpart::rpart(y ~ x)
    purged.model <- purge(unpurged.model)
    test.new.data <- data.frame(x=1:10)
    expect_is(purged.model, 'rpart')
    purge_test_helper(unpurged.model, purged.model, test.new.data)
  }
})

test_that("randomForest purge works correctly", {
  if (requireNamespace('randomForest', quietly=TRUE)) {
    sample.size <- 1000
    x <- rnorm(sample.size)
    y <- x + rnorm(sample.size)
    unpurged.model <- randomForest::randomForest(y ~ x, ntree=10)
    purged.model <- purge(unpurged.model)
    test.new.data <- data.frame(x=1:10)
    expect_is(purged.model, 'randomForest')
    purge_test_helper(unpurged.model, purged.model, test.new.data)
  }
})

test_that("ranger purge works correctly", {
  if (requireNamespace('ranger', quietly=TRUE)) {
    sample.size <- 1000
    x <- rnorm(sample.size)

    # test classification
    y <- as.factor(runif(sample.size) > 0.5)
    unpurged.model <- ranger::ranger(y ~ x, data.frame(x, y),
                                     num.trees=10, write.forest=TRUE)
    purged.model <- purge(unpurged.model)
    test.new.data <- data.frame(x=1:10)
    expect_is(purged.model, 'ranger')
    purge_test_helper(unpurged.model, purged.model, test.new.data,
                      predict.method=function(ranger.model, test.data) {
                        return(predict(ranger.model, test.data)$predictions)
                      })

    # test regression
    y <- rnorm(sample.size)
    unpurged.model <- ranger::ranger(y ~ x, data.frame(x, y),
                                     num.trees=10, write.forest=TRUE)
    purged.model <- purge(unpurged.model)
    expect_is(purged.model, 'ranger')
    purge_test_helper(unpurged.model, purged.model, test.new.data,
                      predict.method=function(ranger.model, test.data) {
                        return(predict(ranger.model, test.data)$predictions)
                      })

    # test survival
    if (requireNamespace('survival', quietly=TRUE)) {
      y.time <- abs(rnorm(sample.size))
      y.status <- ifelse(runif(sample.size) > 0.5, 0, 1)
      unpurged.model <- ranger::ranger(survival::Surv(y.time, y.status) ~ x,
                                       data.frame(x, y),
                                       num.trees=10, write.forest=TRUE)
      purged.model <- purge(unpurged.model)
      expect_is(purged.model, 'ranger')
      purge_test_helper(unpurged.model, purged.model, test.new.data,
                        predict.method=function(ranger.model, test.data) {
                          preds <- predict(ranger.model, test.data)$chf
                          return(preds[, ncol(preds)])
                        })
    }
  }
})

test_that("coxph purge works correctly", {
  if (requireNamespace('survival', quietly=TRUE)) {
    sample.size <- 1000
    x <- rnorm(sample.size)
    y.time <- abs(rnorm(sample.size))
    y.status <- ifelse(runif(sample.size) > 0.5, 0, 1)
    unpurged.model <- survival::coxph(survival::Surv(y.time, y.status) ~ x)
    purged.model <- purge(unpurged.model)
    test.new.data <- data.frame(x=1:10)
    expect_is(purged.model, 'coxph')
    purge_test_helper(unpurged.model, purged.model, test.new.data)
  }
})
