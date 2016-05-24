library(testthat)
library(ff)

context("cut")

test_that("Cut works for ff vectors",{
	x <- 1:10 
	xf <- ff(x)
	expect_identical( cut(x, breaks=3)
	                , cut(xf, breaks=3)[]
				    )
})
