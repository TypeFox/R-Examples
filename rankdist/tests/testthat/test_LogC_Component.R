library(rankdist)
context("LogC_Component")

lst=list(c(0,0),c(0,1),c(1,0),c(2,0),c(1,1),c(2,1))

test_that("correct ouput",{
    expect_equal(exp(LogC_Component(c(1,1))),Reduce(`+`,lapply(lst,function(x)exp(-x[1]-x[2]))))
    expect_equal(exp(LogC_Component(c(1,2))),Reduce(`+`,lapply(lst,function(x)exp(-x[1]-2*x[2]))))
})