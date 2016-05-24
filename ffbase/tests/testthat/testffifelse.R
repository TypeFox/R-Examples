library(testthat)
library(ff)

context("ffifelse")

oldffmaxbytes <- getOption("ffmaxbytes")
options(ffmaxbytes = 20)

test_that("ffifelse works",{
	data(iris)
	ffiris <- as.ffdf(iris)
	
	test.ff <- ffifelse(ffiris$Sepal.Length < 5, TRUE, NA)
	test.ram <- ifelse(iris$Sepal.Length < 5, TRUE, NA)	
	expect_equal(test.ram, test.ff[])
})
test_that("ffifelse works",{
	data(iris)
	ffiris <- as.ffdf(iris)

	test.ff <- ffifelse(ffiris$Sepal.Length < 5, factor("abc"), NA)
	test.ram <- ifelse(iris$Sepal.Length < 5, "abc", NA)	
	expect_equal(test.ram, as.character(test.ff[]))
})
test_that("ffifelse works",{
	data(iris)
	ffiris <- as.ffdf(iris)
	test.ff <- ffifelse(ffiris$Sepal.Length < 5, Sys.Date(), factor("abc"))
	test.ram <- ifelse(iris$Sepal.Length < 5, Sys.Date(), "abc")
	expect_equal(test.ram, as.character(test.ff[]))
})

options(ffmaxbytes = oldffmaxbytes)



