library(rankdist)
context("fai and weight")

fai1 = c(0,0,0,1)
w1 = rep(1,4)
fai2 = rep(0.1,4)
w2 = c(0.4,0.3,0.2,0.1)

test_that("fai to w is correct",{
    expect_equal(paramTow(fai1),w1)
    expect_equal(paramTow(fai2),w2)
})

test_that("w to fai is correct",{
    expect_equal(wToparam(w1),fai1)
    expect_equal(wToparam(w2),fai2)
})