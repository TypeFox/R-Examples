library(metricTester)
context("Ensure metrics, nulls and simulations return correct values")

tree <- geiger::sim.bdtree(b=0.1, d=0, stop="taxa", n=50)
set.seed(0)
sim.abundances <- round(rlnorm(5000, meanlog=2, sdlog=1)) + 1
set.seed(0)
cdm <- simulateComm(tree, richness.vector=10:25, abundances=sim.abundances)
abund <- abundanceVector(cdm)
set.seed(0)
coords <- data.frame(lat=rnorm(n=16, mean=20, sd=1), long=rnorm(16, mean=100, sd=1))
dists <- as.matrix(dist(coords, diag=T, upper=T))
row.names(dists) <- row.names(cdm)
colnames(dists) <- row.names(cdm)
prepped <- prepData(tree=tree, picante.cdm=cdm)
metricResults <- calcMetrics(prepped)
theEnd <- dim(metricResults)[2]
prepped <- prepNulls(tree=tree, picante.cdm=cdm, regional.abundance=abund,
	distances.among=dists)
nullResults <- runNulls(prepped)
prepped <- prepSimulations(tree, arena.length=300, mean.log.individuals=2, 
	length.parameter=5000, sd.parameter=50, max.distance=20, proportion.killed=0.2,
	competition.iterations=3)
simulationResults <- runSimulations(prepped)
noPlots <- dim(cdm)[1]
noSp <- length(tree$tip.label)

test_that("metric result dataframe has correct dimensions",
{
	expect_true(dim(metricResults)[1] == noPlots)
	expect_true(dim(metricResults)[2] == theEnd)
})

test_that("all metric results are real numbers",
{
	#apply is.finite over the dataframe. should return all true. check whether any equal
	#false. if they do, then the sum of falses should be more than 0.
	expect_true(sum(apply(metricResults[,3:theEnd], 2, is.finite)==FALSE)==0)
})

test_that("null results list is same length as number of null models",
{
	expect_true(length(nullResults) == length(defineNulls()))
})

test_that("each null result is a matrix of appropriate dimensions",
{
	expect_true(all(lapply(nullResults, class)=="matrix"))
	expect_true(all(as.data.frame(lapply(nullResults, dim))[1,]==noPlots))
	expect_true(all(as.data.frame(lapply(nullResults, dim))[2,]==noSp))
})

test_that("simulation results list is same length as number of simulations",
{
	expect_true(length(simulationResults) == length(defineSimulations()))
})

test_that("each simulation contains an element called 'arena' with non-zero elements",
{
	expect_true(all(as.data.frame(lapply(seq_along(simulationResults), 
		function(x) dim(simulationResults[[x]]$arena))) > 1))
})

test_that("each simulation contains an element called 'regional abundance'",
{
	expect_true(all(as.data.frame(lapply(seq_along(simulationResults), 
		function(x) length(simulationResults[[x]]$regional.abundance))) > 1))
})
