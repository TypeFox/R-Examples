library(rankdist)
context("RankData.initialize")

ranking = matrix(c(1,2,3,4,4,3,2,1,1,3,2,3,1,2,2,2,2,2,2,1),ncol=4,byrow=TRUE)
# ordering not provided
topq = c(3,2,1)

# supplied nobj, nobs, count correct
f01 = function() new("RankData",ranking=ranking,nobj=4,nobs=4,count=c(2,3,2,2,1))
test0_obj = new("RankData",ranking=ranking,nobj=4,count=c(2,3,2,2,1))
test_that("supplied nobj, nobs, count",{
    expect_error(f01())
    expect_equal(test0_obj@nobj,4)
    expect_equal(test0_obj@nobs,10) 
    expect_equal(test0_obj@count,c(2,3,2,2,1))
})

# default nobj, nobs, count correct
test1_obj = new("RankData",ranking=ranking)
test_that("default for nobj, nobs, count",{
    expect_equal(test1_obj@nobj,4)
    expect_equal(test1_obj@nobs,5) 
    expect_equal(test1_obj@count,rep(1,5))
})

test2_obj = new("RankData",ranking=ranking,topq=topq)
test_that("sequential topq",{
    expect_equal(test2_obj@nobj,4)
    expect_equal(test2_obj@nobs,5) 
    expect_equal(test2_obj@count,rep(1,5))
    expect_equal(test2_obj@topq,c(3,2,1))
    expect_equal(test2_obj@q_ind,c(1,3,4,6))
    expect_equal(test2_obj@subobs,c(2,1,2))
})

ranking = matrix(c(1,2,2,2,2,2,2,1,1,2,3,4,4,3,2,1,1,3,2,3),ncol=4,byrow=TRUE)
topq = c(1,3,2)
test3_obj = new("RankData",ranking=ranking,topq=topq)
f3 = function() new("RankData",ranking=ranking,topq=c(1,4,2))
test_that("non-sequential topq",{
    expect_warning(f3(),"topq value should range between 1 and nobs-1")
    expect_equal(test3_obj@nobj,4)
    expect_equal(test3_obj@nobs,5) 
    expect_equal(test3_obj@count,rep(1,5))
    expect_equal(test3_obj@topq,c(1,3,2))
    expect_equal(test3_obj@q_ind,c(1,3,5,6))
    expect_equal(test3_obj@subobs,c(2,2,1))
})
