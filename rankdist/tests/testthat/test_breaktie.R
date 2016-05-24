library(rankdist)
context("break tie")

obj4 = rbind(1:4, permute::allPerms(4))
p2 = c(3,2,1,3)
p1 = c(1,2,2,2)
count2 = c(rep(1,24), 2)
count1 = c(rep(1,24), 2, 3)
dat2 = new("RankData", ranking=rbind(obj4,p2), count=count2, topq=c(3,2))
res2 = BreakTieEqualProb(dat2)

dat2_p = new("RankData", ranking=matrix(p2,nrow=1), count=2, topq=2, nobj=4)
res2_p = BreakTieEqualProb(dat2_p)

dat1 = new("RankData", ranking=rbind(obj4, p2, p1), count=count1, topq=c(3,2,1))
res1 = BreakTieEqualProb(dat1)
which(res1@count == 3)


test_that("break tie is correct",{
    expect_equal(which(res2@count==2), c(15,21))
    expect_equal(res2_p@ranking,rbind(c(3,2,1,4), c(4,2,1,3)))
    expect_equal(which(res1@count == 1.5), 1:6)
    expect_equal(which(res1@count == 2), c(15,21))
})
