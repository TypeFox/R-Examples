context("odd.prob")

test_that("odd.prob outputs correctly", {
    odds <- 1
    expect_true(odd.prob(odds) == 0.5)
}) 
