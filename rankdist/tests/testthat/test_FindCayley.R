library(rankdist)
context("FindCayley")

mat = rbind(c(4,2,1,3), c(1,2,3,4), c(4,3,2,1))
mat2 = rbind(c(2,5,3,1,4,5))
pi0 = c(1,2,3,4)
pi02 = 1:6

test_that("correct ouput",{
    expect_equal(FindCayley(mat, pi0), c(2,0,2))
    expect_equal(FindCayley(mat, 4:1), c(2,2,0))
    expect_equal(FindCayley(mat2, pi02), 3)
})
