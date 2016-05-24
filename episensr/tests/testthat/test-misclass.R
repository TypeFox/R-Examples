context("check misclassification bias")

test_that("correct number of arguments for bias parameters", {
    expect_that(misclassification(matrix(c(215, 1449, 668, 4296),
                                         nrow = 2, byrow = TRUE),
                                  implement = "exposure",
                                  bias = c(.99, .99)),
                throws_error())
})

test_that("bias parameters between 0 and 1", {
    expect_that(misclassification(matrix(c(215, 1449, 668, 4296),
                                         nrow = 2, byrow = TRUE),
                                  implement = "exposure",
                                  bias = c(-1, .78, .99, 2)),
                throws_error())
})

test_that("Observed measures are correct for exposure misclassification", {
    model <- misclassification(matrix(c(215, 1449, 668, 4296),
                                         nrow = 2, byrow = TRUE),
                                  implement = "exposure",
                                  bias = c(.78, .78, .99, .99), print = FALSE)
    expect_equal(model$obs.measures[1, 1], 0.9654, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 0.8524, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 1.0933, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 1], 0.9542, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 2], 0.8093, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 3], 1.1252, tolerance = 1e-4, scale = 1)
})

test_that("Adjusted measures are correct for exposure misclassification", {
    model <- misclassification(matrix(c(215, 1449, 668, 4296),
                                         nrow = 2, byrow = TRUE),
                                  implement = "exposure",
                                  bias = c(.78, .78, .99, .99), print = FALSE)
    expect_equal(model$adj.measures[1], 0.9614, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2], 0.9491, tolerance = 1e-4, scale = 1)
})

test_that("Observed measures are correct for outcome misclassification", {
    model <- misclassification(matrix(c(4558, 3428, 46305, 46085),
                                         nrow = 2, byrow = TRUE),
                                  implement = "outcome",
                                  bias = c(.53, .53, .99, .99), print = FALSE)
    expect_equal(model$obs.measures[1, 1], 1.2944, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 1.2404, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 1.3506, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 1], 1.3233, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 2], 1.2636, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 3], 1.3858, tolerance = 1e-4, scale = 1)
})

test_that("Adjusted measures are correct for outcome misclassification", {
    model <- misclassification(matrix(c(4558, 3428, 46305, 46085),
                                         nrow = 2, byrow = TRUE),
                                  implement = "outcome",
                                  bias = c(.53, .53, .99, .99), print = FALSE)
    expect_equal(model$adj.measures[1], 1.3440, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2], 1.4062, tolerance = 1e-4, scale = 1)
})
