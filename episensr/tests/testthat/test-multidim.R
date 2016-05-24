context("check multidimensional bias")

test_that("Observed measures are correct for exposure misclassification", {
    model <- multidimBias(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
             type = "exposure",
             se = c(1, 1, 1, .9, .9, .9, .8, .8, .8),
             sp = c(1, .9, .8, 1, .9, .8, 1, .9, .8), print = FALSE)
    expect_equal(model$obs.measures[1, 1], 1.6470, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 1.1824, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 2.2941, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 1], 1.7603, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 2], 1.2025, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 3], 2.5769, tolerance = 1e-4, scale = 1)
})

test_that("Adjusted measures are correct for exposure misclassification", {
    model <- multidimBias(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
             type = "exposure",
             se = c(1, 1, 1, .9, .9, .9, .8, .8, .8),
             sp = c(1, .9, .8, 1, .9, .8, 1, .9, .8), print = FALSE)
    expect_equal(model$adj.measures[[1]][1, 1], 1.6470, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[[2]][1, 1], 1.7603, tolerance = 1e-4, scale = 1)
})

test_that("Observed measures are correct for outcome misclassification", {
    model <- multidimBias(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
             type = "exposure",
             se = c(1, 1, 1, .9, .9, .9, .8, .8, .8),
             sp = c(1, .9, .8, 1, .9, .8, 1, .9, .8), print = FALSE)
    expect_equal(model$obs.measures[1, 1], 1.6470, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 1.1824, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 2.2941, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 1], 1.7603, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 2], 1.2025, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 3], 2.5769, tolerance = 1e-4, scale = 1)
})

test_that("Adjusted measures are correct for outcome misclassification", {
    model <- multidimBias(matrix(c(45, 94, 257, 945), nrow = 2, byrow = TRUE),
             type = "outcome",
             se = c(1, 1, 1, .9, .9, .9, .8, .8, .8),
             sp = c(1, .9, .8, 1, .9, .8, 1, .9, .8), print = FALSE)
    expect_equal(model$adj.measures[[1]][1, 1], 1.6470, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[[2]][1, 1], 1.7603, tolerance = 1e-4, scale = 1)
})

test_that("Observed measures are correct for confounder misclassification", {
    model <- multidimBias(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
             type = "confounder",
             bias = list(seq(.72, .92, by = .02),
                 seq(.01, .11, by = .01), seq(.13, 1.13, by = .1)), print = FALSE)
    expect_equal(model$obs.measures[1, 1], 0.3479, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 0.2757, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 0.4390, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 1], 0.2180, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 2], 0.1519, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 3], 0.3128, tolerance = 1e-4, scale = 1)
})

test_that("Adjusted measures are correct for confounder misclassification", {
    model <- multidimBias(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
             type = "confounder",
             bias = list(seq(.72, .92, by = .02),
                 seq(.01, .11, by = .01), seq(.13, 1.13, by = .1)), print = FALSE)
    expect_equal(model$adj.measures[1, 1], 0.9232, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 0.7748, tolerance = 1e-4, scale = 1)
})

test_that("Observed measures are correct for selection misclassification", {
    model <- multidimBias(matrix(c(136, 107, 297, 165), nrow = 2, byrow = TRUE),
             type = "selection",
             OR.sel = seq(1.5, 6.5, by = .5), print = FALSE)
    expect_equal(model$obs.measures[1, 1], 0.7984, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 0.6518, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 0.9780, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 1], 0.7061, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 2], 0.5144, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 3], 0.9693, tolerance = 1e-4, scale = 1)
})

test_that("Adjusted measures are correct for selection misclassification", {
    model <- multidimBias(matrix(c(136, 107, 297, 165), nrow = 2, byrow = TRUE),
             type = "selection",
             OR.sel = seq(1.5, 6.5, by = .5), print = FALSE)
    expect_equal(model$adj.measures[1, 2], 0.4708, tolerance = 1e-1, scale = 1,
                 check.attributes = FALSE)
    expect_equal(model$adj.measures[2, 2], 0.3531, tolerance = 1e-1, scale = 1,
                 check.attributes = FALSE)
})
