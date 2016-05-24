library(mlsjunkgen)
context("mlsjunkgend")

dat <- data.frame(0.95516)
names(dat) <- "RN"
dat2 <- data.frame(c(0.95516, 0.66908, 0.21235, 0.34488, 0.11995))
names(dat2) <- "RN"
dat3 <- data.frame(c(0.887, 0.762, 0.214, 0.197, 0.375, 0.168))
names(dat3) <- "RN"


test_that("mlsjunkgend returns the correct data frame", {
    expect_equal(mlsjunkgend(w = 1, x = 2, y = 3, z = 4), dat)
    expect_equal(mlsjunkgend(n = 5, w = 1, x = 2, y = 3, z = 4), dat2)
    expect_equal(mlsjunkgend(n = 6, w = 2, x = 3, y = 4, z = 5, round = 3), 
                 dat3)
    expect_error(mlsjunkgend(n = "X", w = 1, x = 2, y = 3, z = 4), "Invalid input.  Please ensure n is numeric.")
    expect_error(mlsjunkgend(n = 5, w = "X", x = 2, y = 3, z = 4), "Invalid input.  Please ensure all seeds are numeric.")
})