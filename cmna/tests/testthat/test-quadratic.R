context("quadratic")

## Someone of these should return complex, but they don't because R does
## not support complex natively.
expect_equal(quadratic(1, 0, -1), c(-1, 1))
expect_equal(quadratic(4, -4, 1), c(0.5, 0.5))

expect_equal(quadratic2(1, 0, -1), c(-1, 1))
expect_equal(quadratic2(4, -4, 1), c(0.5, 0.5))
