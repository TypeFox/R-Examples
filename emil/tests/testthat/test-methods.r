context("Methods")

test_that("pamr", {
    reset_notification()
    expect_message(fit <- fit_pamr(iris[-5], iris$Species,
            cv=resample("crossvalidation", iris$Species, nfold=5, nrepeat=1)),
        "Use.*pre_pamr.*pre-processing")
    expect_that(fit, is_a("list"))
})

test_that("glmnet", {
    reset_notification()
    expect_message(fit <- fit_glmnet(iris[-5], iris$Species),
                   ".*data set.*matrix.*form.*")
    model <- fit_glmnet(as.matrix(iris[-5]), iris$Species)
    expect_that(fit, is_a("list"))
    expect_that(fit$glmnet.fit, is_a("glmnet"))
})
