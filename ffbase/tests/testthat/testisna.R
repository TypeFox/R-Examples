library(testthat)
library(ff)

context("isna")

test_that("isna works",{
	size <- 100
	x <- rnorm(size)
	x[sample(1:size, round(size/2))] <- NA
	test.ram <- is.na(x)
	test.ff <- is.na.ff(ff(x), by=2)
	expect_equal( test.ram, test.ff[])
	
	x <- rnorm(size)
	x.ff <- as.ff(x)
	idx <- sample(1:size, round(size/2))
	is.na(x.ff) <- idx
	x[idx] <- NA
	test.ram <- is.na(x)
	test.ff <- is.na(x.ff)
	expect_equal( test.ram, test.ff[])

	x <- rnorm(size)
	x.ff <- as.ff(x)
	idx <- sample(1:size, round(size/2))
	x[idx] <- NA
	idx <- as.ff(idx)
	is.na(x.ff) <- idx
	test.ram <- is.na(x)
	test.ff <- is.na(x.ff)
	expect_equal( test.ram, test.ff[])
})
