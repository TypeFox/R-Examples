library(metricTester)
context("Complicated community simulations and component functions for final sims")

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
rawResults <- metricsNnulls(tree=tree, picante.cdm=cdm, regional.abundance=abund,
	distances.among=dists, randomizations=3)
results <- reduceRandomizations(rawResults)
observed <- observedMetrics(tree, cdm)
test1 <- errorChecker(observed, results, "richness")
test2 <- errorChecker(observed, results, "both")

test_that("metricsNnulls output a list of randomizations of appropriate length",
{
	expect_true(length(rawResults) == 3)
	#this confirms that each randomization makes a list of data frames of length equal
	#to number of null models
	expect_true(all(lapply(rawResults, length)==length(defineNulls())))
	#this confirms (just for the first randomization) that all metrics are calculated for
	#each null model
	expect_true(all(as.data.frame(lapply(rawResults[[1]], dim))[2,]
		==length(defineMetrics())+1))
})

test_that("reduceRandomizations boils raw results down into appropriate dimensions",
{
	expect_true(length(results) == length(defineNulls()))
	expect_true(all(as.data.frame(lapply(results, dim))[2,]==length(defineMetrics())+1))
})

test_that("errorChecker gives appropriate outputs if concat by equals 'richness'",
{
	expect_true(length(test1) == 2)
	expect_true(all(names(test1$ses) == names(defineNulls())))
})

