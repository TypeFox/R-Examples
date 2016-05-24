library(testthat)

context("droplevels")

test_that("droplevels.ff works",{
    x <- factor(c("c","b"), levels=c("a","b","c"))
	fx <- ff(x)
	expect_identical(droplevels(x), (droplevels(fx))[])
	expect_identical(droplevels(x), (fx <- droplevels(fx, inplace=TRUE))[])
})


test_that("droplevels.ffdf works",{
    x <- factor(c("c","b"), levels=c("a","b","c"))
 
	dat <- data.frame(x=x, y=x, z=x)
	fdat <- as.ffdf(dat)
	
	expect_identical(droplevels(dat), (droplevels(fdat))[,]) 
	expect_identical(droplevels(dat,except="x"), (droplevels(fdat, except="x"))[,]) 
})

test_that("droplevels.ffdf keeps names",{
  x <- data.frame(a = 1:10, b = 1:10, c = factor(LETTERS[1:10], levels=LETTERS))
  x <- as.ffdf(x)
  
  names(x) <- c("x","y","z")
  expect_identical(names(droplevels(x)), c("x","y","z"))
})
