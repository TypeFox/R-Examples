context("logit.prob")

test_that("logit.prob outputs correctly", {
    odds <- 1
    logit <- log(odds)
    expect_true(logit.prob(logit) == 0.5)
}) 
