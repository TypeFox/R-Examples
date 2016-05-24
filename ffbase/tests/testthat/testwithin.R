library(testthat)
library(ff)

context("within")

test_that("Within works",{
	dat <- data.frame(x=1:10, y=10:1) 
	ffdat <- as.ffdf(dat)
   
  dat <- within(dat, z<-x+y)
  ffdat <- within(ffdat, z<-x+y)
      
	expect_equivalent( dat
	                 , ffdat[,]
				     )
})