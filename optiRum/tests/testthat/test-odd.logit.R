context("odd.logit")

test_that("odd.logit outputs correctly", {
    odds <- 1
    expect_true(odd.logit(odds) == 0)
    check <- log(0.5)
    odds <- c(1, 0.5, 2)
    expect_true(identical(odd.logit(odds), c(0, check, -check)))
}) 
