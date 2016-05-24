context("horner")

## If you're curious, it's from the Addams Family movie
beta <- c(2, 10, 11)
expect_equal(horner(5, beta), 327)
expect_equal(horner(87, beta), 84131)
expect_equal(horner(-10, beta), 1002)
expect_equal(horner(0, beta), 2)

beta <- c(-1, 0, 1)
expect_equal(horner(c(1, 2, 3), beta), c(0, 3, 8))
