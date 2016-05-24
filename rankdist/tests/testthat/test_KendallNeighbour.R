library(rankdist)
context("Kendall Neighbour")

rank = 1:5


test_that("Output type is correct",{
    expect_is(KendallNeighbour(rank),"matrix")
    expect_equal(nrow(KendallNeighbour(rank)),4)
    expect_equal(ncol(KendallNeighbour(rank)),5)
})


test_that("Output value is correct",{
    expect_equal_to_reference(file="KendallNeighbour.rds",KendallNeighbour(rank))
})