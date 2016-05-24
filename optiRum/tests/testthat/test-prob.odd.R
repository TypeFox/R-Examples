context("prob.odd")

test_that("prob.odd outputs correctly", {
    probs <- 0.5
    expect_true(prob.odd(probs) == 1)
    probs <- c(0.5, 0.25)
    expect_true(identical(prob.odd(probs), c(1, 1/3)))
}) 
