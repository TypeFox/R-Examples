context("prob.logit")

test_that("prob.logit outputs correctly", {
    probs <- 0.5
    expect_true(prob.logit(probs) == 0)
}) 
