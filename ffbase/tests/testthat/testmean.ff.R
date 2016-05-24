library(testthat)
library(ff)


test_that("Mean works",{
	x <- runif(100)
	fx <- ff(x)
	expect_equal(mean(x), mean(fx))
})