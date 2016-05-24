context("check bias due to confounders with effect measure modification")

test_that("correct number of arguments for prevalence of the confounder", {
    expect_that(confounders.emm(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                            implement = "RR",
                            p = c(.8),
                            RR.cd = c(.4, .7)),
                throws_error())
})

test_that("confounder prevalences between 0 and 1", {
    expect_that(confounders.emm(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                            implement = "RR",
                            p = c(-1, 2, .2, .2),
                            RR.cd = c(.4, .7)),
                throws_error())
})

test_that("RR is correct", {
    model <- confounders.emm(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                         implement = "RR",
                         p = c(.8, .05),
                         RR.cd = c(.4, .7), print = FALSE)
    expect_equal(model$obs.measures[1, 1], 0.3479, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 0.2757, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 0.4390, tolerance = 1e-4, scale = 1)
})

test_that("Adjusted RR are correct", {
    model <- confounders.emm(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                         implement = "RR",
                         p = c(.8, .05),
                         RR.cd = c(.4, .7), print = FALSE)
    expect_equal(model$adj.measures[1, 1], 0.4509, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 0.7716, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 0.6370, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.5462, tolerance = 1e-4, scale = 1)
})

test_that("OR is correct", {
    model <- confounders.emm(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                         implement = "OR",
                         p = c(.8, .05),
                         OR.cd = c(.4, .7), print = FALSE)
    expect_equal(model$obs.measures[1, 1], 0.2180, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], 0.1519, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], 0.3128, tolerance = 1e-4, scale = 1)
})

test_that("Adjusted OR are correct", {
    model <- confounders.emm(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                         implement = "OR",
                         p = c(.8, .05),
                         OR.cd = c(.4, .7), print = FALSE)
    expect_equal(model$adj.measures[1, 1], 0.2825, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], 0.7716, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 1], 0.3969, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[2, 2], 0.5492, tolerance = 1e-4, scale = 1)
})

test_that("RD is correct", {
    model <- confounders.emm(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                         implement = "RD",
                         p = c(.8, .05),
                         RD.cd = c(-.6, -.3), print = FALSE)
    expect_equal(model$obs.measures[1, 1], -0.3114, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 2], -0.3903, tolerance = 1e-4, scale = 1)
    expect_equal(model$obs.measures[1, 3], -0.2325, tolerance = 1e-4, scale = 1)
})

test_that("Adjusted RD are correct", {
    model <- confounders.emm(matrix(c(105, 85, 527, 93), nrow = 2, byrow = TRUE),
                         implement = "RD",
                         p = c(.8, .05),
                         RD.cd = c(-.6, -.3), print = FALSE)
    expect_equal(model$adj.measures[1, 1], 0.1213, tolerance = 1e-4, scale = 1)
    expect_equal(model$adj.measures[1, 2], -0.4326, tolerance = 1e-4, scale = 1)
})
