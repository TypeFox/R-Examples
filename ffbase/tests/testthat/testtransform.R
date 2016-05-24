library(testthat)
library(ff)

context("transform")

test_that("transform works",{
	dat <- data.frame(x=1:10, y=10:1) 
	ffdat <- as.ffdf(dat)
   
    dat <- transform(dat, z=x+y)
    ffdat <- transform(ffdat, z=x+y)
      
	expect_equivalent( dat
	                 , ffdat[,]
				     )
})