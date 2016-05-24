context("scaledScore")

test_that("scaledscore tranforms logits correctly", {
    odds <- 1
    logit <- log(odds)
    expect_true(scaledScore(logit) == 300)
    odds <- 0.5
    logit <- log(odds)
    expect_true(scaledScore(logit) == 250)
    odds <- 0.25
    logit <- log(odds)
    expect_true(scaledScore(logit) == 200)
    odds <- 2
    logit <- log(odds)
    expect_true(scaledScore(logit) == 350)
    odds <- 4
    logit <- log(odds)
    expect_true(scaledScore(logit) == 400)
    odds <- 0.33
    logit <- log(odds)
    expect_true(scaledScore(logit) == 220)
})

test_that("scaledscore tranforms logits correctly", {
    odds <- c(1, 0.5, 0.25, 2, 4)
    logit <- log(odds)
    expect_true(identical(scaledScore(logit), c(300, 250, 200, 350, 400)))
}) 
