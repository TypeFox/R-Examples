library(rankdist)
context("Cayley Neighbour")

rank = 1:4

test_that("Output type is correct",{
    expect_is(CayleyNeighbour(rank),"matrix")
    expect_equal(nrow(CayleyNeighbour(rank)),6)
    expect_equal(ncol(CayleyNeighbour(rank)),4)
})

test_that("Output value is correct",{
    expect_equal_to_reference(file="CayleyNeighbour.rds",CayleyNeighbour(rank))
})