context("check limits on confounding")

test_that("Limits are correct", {
    model <- confounders.limit(OR = 1.65, crude.RR = 1.5, print = FALSE)
    expect_equal(model$conf.limits[1], 0.9091, tolerance = 1e-4, scale = 1)
    expect_equal(model$conf.limits[2], 1.5000, tolerance = 1e-4, scale = 1)
})
