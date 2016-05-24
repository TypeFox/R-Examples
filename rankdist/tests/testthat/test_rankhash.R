library(rankdist)
context("RankHash")

rank = sample(1:52,52)
rankmat = replicate(10,sample(1:52,52), simplify = "array")

test_that("Rank Hash",{
    expect_equal(HashtoRank(RanktoHash(rank)),rank)
    expect_equal(HashtoRank(RanktoHash(HashtoRank(RanktoHash(rank)))),rank)
    expect_equal(HashtoRank(RanktoHash(rankmat)),rankmat)
    expect_equal(HashtoRank(RanktoHash(HashtoRank(RanktoHash(rankmat)))),rankmat)
})