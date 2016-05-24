library(rankdist)
context("FindV")

mat=matrix(ncol=4,c(3,4,1,2,2,1,4,3,1,2,3,4),byrow=TRUE)
mat2 = mat
pi0 = c(3,4,1,2)
pi02 = pi0
resmat = matrix(ncol=3,c(0,0,0,3,2,1,2,2,0),byrow=TRUE)
test_that("correct ouput",{
    expect_equal(FindV(mat,pi0),resmat)
    expect_equal(mat,mat2) # no side effect
    expect_equal(pi0,pi02)
})

