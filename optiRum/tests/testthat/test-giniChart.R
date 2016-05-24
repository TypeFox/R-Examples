context("giniChart")
library("ggplot2")
test_that("giniChart correctly produce a chart, numeric outcome", {
    rm(list = ls())
    sampledata <- data.frame(val = rnorm(100), outcome = rbinom(100, 1, 0.8))
    check1 <- ggplot_build(giniChart(sampledata$val, sampledata$outcome))
    expect_that(check1, is_a("list"))
})

test_that("giniChart correctly produce a chart, factor outcome", {
    rm(list = ls())
    sampledata <- data.frame(val = rnorm(100), outcome = factor(rbinom(100, 1, 0.8)))
    check1 <- ggplot_build(giniChart(sampledata$val, sampledata$outcome))
    expect_that(check1, is_a("list"))
})

test_that("giniChart errors given incorrect input to pred", {
    rm(list = ls())
    sampledata <- data.frame(val = rnorm(100), outcome = factor(rbinom(100, 1, 0.8)))
    expect_error(giniChart(sampledata$outcome, sampledata$val))
})

test_that("giniChart errors given continuous input to act", {
    rm(list = ls())
    sampledata <- data.frame(val = rnorm(100), outcome = factor(rbinom(100, 1, 0.8)))
    expect_error(giniChart(sampledata$val, sampledata$val))
})

test_that("giniChart warns given character input to act", {
    rm(list = ls())
    sampledata <- data.frame(val = rnorm(26), outcome = letters, stringsAsFactors = FALSE)
    expect_error(giniCoef(sampledata$val, sampledata$outcome))
}) 
