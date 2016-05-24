library(testthat)
library(ff)

context("duplicated")

test_that("duplicated.ff works",{
  data(iris)
	irisdouble <- rbind(iris, iris)
	irisdouble <- irisdouble[sample(x=1:nrow(irisdouble), size=nrow(irisdouble), replace = FALSE), ]
	rownames(irisdouble) <- NULL
	
	ffiris <- as.ffdf(irisdouble)
	test.ff <- duplicated(ffdfsort(ffiris), by=10)
	test.ram <- duplicated(irisdouble[order(irisdouble[, 1], irisdouble[, 2], irisdouble[, 3], irisdouble[, 4], irisdouble[, 5]), ])
	expect_equal( test.ram, test.ff[])
	
	test.ff <- duplicated(ffsort(ffiris$Sepal.Length), by=10)
	test.ram <- duplicated(irisdouble$Sepal.Length[order(irisdouble$Sepal.Length)])
	expect_equal( test.ram, test.ff[])
	  
})

