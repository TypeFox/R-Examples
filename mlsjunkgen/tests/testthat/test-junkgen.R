library(mlsjunkgen)
context("junkgen")

test_that("junkgen returns the correct psuedo-random number", {
    expect_equal(round(junkgen(w = 1, x = 2, y = 3, z = 4), 5), 0.95516)
    expect_equal(round(junkgen(w = 2, x = 3, y = 4, z = 5), 5), 0.88696)
    expect_equal(round(junkgen(w = 0.886956, x = 0.886956, y = 0.9551644, 
                               z = 0.9551644), 5), 0.32323)
    expect_error(junkgen(w = "X", x = 2, y = 3, z = 4), "Invalid input.  Please ensure all seeds are numeric.")
})