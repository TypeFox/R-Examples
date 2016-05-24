context("giniCoef")

test_that("giniCoef correctly produces a value", {
    rm(list = ls())
    sampledata <- data.frame(val = rnorm(100), outcome = rbinom(100, 1, 0.8))
    check1 <- giniCoef(sampledata$val, sampledata$outcome)
    expect_true(abs(check1) <= 1)
})

test_that("giniCoef correctly produce a value, factor outcome", {
    rm(list = ls())
    sampledata <- data.frame(val = rnorm(100), outcome = factor(rbinom(100, 1, 0.8)))
    check1 <- giniCoef(sampledata$val, sampledata$outcome)
    expect_true(abs(check1) <= 1)
})


test_that("giniCoef errors given incorrect input to pred", {
    rm(list = ls())
    sampledata <- data.frame(val = rnorm(100), outcome = factor(rbinom(100, 1, 0.8)))
    expect_error(giniCoef(sampledata$outcome, sampledata$val))
})

test_that("giniCoef warns given continuous input to act", {
    rm(list = ls())
    sampledata <- data.frame(val = rnorm(100), outcome = factor(rbinom(100, 1, 0.8)))
    expect_error(giniCoef(sampledata$val, sampledata$sampledata))
})

test_that("giniCoef errors given character input to act", {
    rm(list = ls())
    sampledata <- data.frame(val = rnorm(26), outcome = letters, stringsAsFactors = FALSE)
    expect_error(giniCoef(sampledata$val, sampledata$outcome))
}) 
