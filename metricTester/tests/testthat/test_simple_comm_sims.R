library(metricTester)
context("Basic community simulations and utility functions")

tree <- geiger::sim.bdtree(b=0.1, d=0, stop="taxa", n=50)
set.seed(0)
sim.abundances <- round(rlnorm(5000, meanlog=2, sdlog=1)) + 1
set.seed(0)
cdm <- simulateComm(tree, richness.vector=10:25, abundances=sim.abundances)
abund <- abundanceVector(cdm)

test_that("cdm has correct dimensions",
{
	expect_true(dim(cdm)[1] == 16)
	expect_true(dim(cdm)[2] == 50)
})

test_that("cdm is class data frame",
{
	expect_is(cdm, "data.frame")
})

test_that("cdm properly filled",
{
	expect_true(sum(cdm, na.rm=TRUE) != 0)
})

test_that("row and column names are correct",
{
	expect_true(sum(row.names(cdm) == paste("plot",1:dim(cdm)[1], sep=""))==16)
	expect_true(sum(names(cdm) == tree$tip.label)==50)
})

test_that("abundance vector is vector of class character",
{
	expect_is(abund, "character")
})

test_that("abundance vector has length > 1",
{
	expect_true(length(abund) > 1)
})
