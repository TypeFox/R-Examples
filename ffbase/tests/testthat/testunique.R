library(testthat)
library(ff)

context("unique")

test_that("unique works",{
  data(iris)
	irisdouble <- rbind(iris, iris)
	irisdouble[sample(1:nrow(irisdouble), round(nrow(irisdouble)/10)), ] <- NA
	ffiris <- as.ffdf(irisdouble)
	test.ff <- unique(ffiris$Sepal.Length, by=10)
	test.ram <- unique(ffiris$Sepal.Length[])
	
	expect_true(!sum(!test.ff[] %in% test.ram))
	expect_true(!sum(!test.ram %in% test.ff[]))

	test.ff <- unique(ffiris[1:3], by=10)
	test.ram <- unique(ffiris[, 1:3][,])
 
  expect_true(!sum(!apply(test.ff[,], MARGIN=1, FUN=function(x) paste(x, collapse=",")) %in% apply(test.ram, MARGIN=1, FUN=function(x) paste(x, collapse=","))))
	expect_true(!sum(!apply(test.ram, MARGIN=1, FUN=function(x) paste(x, collapse=",")) %in% apply(test.ff[,], MARGIN=1, FUN=function(x) paste(x, collapse=","))))
  
})

