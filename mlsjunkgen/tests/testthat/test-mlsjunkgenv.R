library(mlsjunkgen)
context("mlsjunkgenv")

test_that("mlsjunkgenv returns the correct vector", {
    expect_equal(mlsjunkgenv(w = 1, x = 2, y = 3, z = 4), 0.95516)
    expect_equal(mlsjunkgenv(n = 5, w = 1, x = 2, y = 3, z = 4), 
                 c(0.95516, 0.66908, 0.21235, 0.34488, 0.11995))
    expect_equal(mlsjunkgenv(n = 6, w = 2, x = 3, y = 4, z = 5, round = 3), 
                 c(0.887, 0.762, 0.214, 0.197, 0.375, 0.168))
    expect_error(mlsjunkgenv(n = "X", w = 1, x = 2, y = 3, z = 4), "Invalid input.  Please ensure n is numeric.")
    expect_error(mlsjunkgenv(n = 5, w = "X", x = 2, y = 3, z = 4), "Invalid input.  Please ensure all seeds are numeric.")
})