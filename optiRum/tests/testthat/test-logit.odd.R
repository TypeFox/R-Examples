context("logit.odd")

test_that("logit.odd outputs correctly", {
    odds <- 1
    logit <- log(odds)
    expect_true(logit.odd(logit) == odds)
    odds <- c(1, 0.5, 0.25, 2, 4)
    logit <- log(odds)
    expect_true(identical(logit.odd(logit), odds))
}) 
