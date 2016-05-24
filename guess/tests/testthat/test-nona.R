context("Test nona")

test_that("nona works correctly", {
	x <- c(NA, 1, 0)
	y <- nona(x)
	expect_that(y, equals(c(0,1,0)))
})