library(testthat)
library(ff)

context("ffmatch")

test_that("Matching based on ffmatch works for ff_vector",{
	x <- as.factor(c(LETTERS, NA))
	y <- as.factor(c("C","B","Z","X","HMM","Nothing",NA))
	test.ram <- match(x, y)
	test.ff <- ffmatch(as.ff(x), as.ff(y), by=2)
	expect_identical( test.ram
	                , test.ff[]
				    )
})

test_that("Matching based on ffmatch works for ffdf", {
  data(iris)
  iris <- unique(iris)
	ffiris <- as.ffdf(iris)
	ffirissubset <- as.ffdf(iris[c(1:10, nrow(iris)), ])
	test.ff <- ffdfmatch(ffiris, ffirissubset, by=2)
	test.ram <- match(apply(ffiris[,], MARGIN=1, FUN=function(x) paste(x, collapse="")), apply(ffirissubset[,], MARGIN=1, FUN=function(x) paste(x, collapse="")))
	
	expect_identical( test.ram
	                , test.ff[]
				    )
})

test_that("The %in% operator based on ffmatch works for ff_vector", {
	x <- as.factor(c(LETTERS, NA))
	y <- as.factor(c("C","B","Z","X","HMM","Nothing",NA))
 	test.ff <- as.ff(x) %in% as.ff(y)
	test.ram <- x %in% y

  expect_equivalent( test.ram
                   , test.ff[]
                   )  
})

test_that("The %in% operator based on ffmatch works for ffdf", {
  data(iris)
  iris <- unique(iris)
	ffiris <- as.ffdf(iris)
	ffirissubset <- as.ffdf(iris[c(1:10, nrow(iris)), ])
	test.ff <- ffiris %in% ffirissubset
	test.ram <- apply(ffiris[,], MARGIN=1, FUN=function(x) paste(x, collapse="")) %in% apply(ffirissubset[,], MARGIN=1, FUN=function(x) paste(x, collapse=""))
	
  expect_equivalent( test.ram
                   , test.ff[]
                   )   
})
