context("check selection bias")

test_that("correct number of arguments for selection probabilities", {
    expect_that(selection(matrix(c(136, 107, 297, 165),
                                 nrow = 2, byrow = TRUE),
                          selprob = c(.94, .85, .64)),
                throws_error())
})

test_that("selection probabilities between 0 and 1", {
    expect_that(selection(matrix(c(136, 107, 297, 165),
                                 nrow = 2, byrow = TRUE),
                          selprob = c(.94, -1, .64, 3)),
                throws_error())
})

test_that("RR is correct", {
    model <- selection(matrix(c(136, 107, 297, 165),
                                 nrow = 2, byrow = TRUE),
                          selprob = c(.94, .85, .64, .25), print = FALSE)
    expect_equal(model$obs.measures[1, 1], 0.7984, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 0.6518, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 0.9780, tolerance = 1e-4, scale = 1)
})

test_that("OR is correct", {
    model <- selection(matrix(c(136, 107, 297, 165),
                                 nrow = 2, byrow = TRUE),
                          selprob = c(.94, .85, .64, .25), print = FALSE)
    expect_equal(model$obs.measures[2, 1], 0.7061, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 2], 0.5144, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[2, 3], 0.9693, tolerance = 1e-4, scale = 1)
})

test_that("Adjusted RR is correct", {
    model <- selection(matrix(c(136, 107, 297, 165),
                                 nrow = 2, byrow = TRUE),
                          selprob = c(.94, .85, .64, .25), print = FALSE)
    expect_equal(model$adj.measures[1], 1.4838, tolerance = 1e-4, scale = 1)
})

test_that("Adjusted OR is correct", {
    model <- selection(matrix(c(136, 107, 297, 165),
                                 nrow = 2, byrow = TRUE),
                          selprob = c(.94, .85, .64, .25), print = FALSE)
    expect_equal(model$adj.measures[2], 1.6346, tolerance = 1e-4, scale = 1)
})
