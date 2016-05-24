library(mlsjunkgen)
context("mlsjunkgenm")

mat <- matrix(0.95516)
mat2 <- matrix(c(0.95516, 0.66908, 0.21235, 0.34488, 0.11995, 0.56398, 0.59235, 
                 0.11432, 0.33525, 0.70271, 0.41810, 0.31337, 0.91985, 0.37872, 
                 0.28042), nrow = 5, ncol = 3)

test_that("mlsjunkgenm returns the correct matrix", {
    expect_equal(mlsjunkgenm(w = 1, x = 2, y = 3, z = 4), mat)
    expect_equal(mlsjunkgenm(nrow = 5, ncol = 3, w = 1, x = 2, y = 3, z = 4), mat2)
    expect_error(mlsjunkgenm(nrow = "X", ncol = 1, w = 1, x = 2, y = 3, z = 4), "Invalid input.  Please ensure nrow and ncol are numeric.")
    expect_error(mlsjunkgenm(w = "X", x = 2, y = 3, z = 4), "Invalid input.  Please ensure all seeds are numeric.")
})