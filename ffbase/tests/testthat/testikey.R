library(testthat)
library(ff)

context("ikey")

test_that("ikey works",{
  
  opts <- options(ffmaxbytes=20)
  on.exit(options(opts))
  
	data(iris)
	ffiris <- iris[dforder(iris[c("Petal.Width","Species")][,]),]
	rownames(ffiris) <- NULL
	ffiris <- as.ffdf(ffiris)
	
	test.ram <- ffiris[c("Petal.Width","Species")][,]
	test.ram <- cumsum(!duplicated(test.ram))
	test.ff <- ikey(ffiris[c("Petal.Width","Species")])
	expect_equal( test.ram
	            , test.ff[]
				      )
})
