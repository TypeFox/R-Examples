library(metricTester)
context("Functions related to sampling plots from arenas")

tree <- geiger::sim.bdtree(b=0.1, d=0, stop="taxa", n=50)
temp <- evolveTraits(tree)
prepped <- prepSimulations(tree, arena.length=300, mean.log.individuals=4, 
	length.parameter=5000, sd.parameter=50, max.distance=20, proportion.killed=0.2,
	competition.iterations=3)
singleArena <- filteringArena(prepped)
bounds <- plotPlacer(no.plots=20, arena.length=300, plot.length=sqrt(1000))
cdm <- plotContents(singleArena$arena, bounds)

#make a simple for loop to run through the bounds and see if any are overlapping.
#make an empty vector to save errors into
error <- c()
for(i in 1:dim(bounds$plot.bounds)[1])
{
	for(j in 1:dim(bounds$plot.bounds)[1])
	{
		#if X1 is bigger than another X1 and less than the corresponding X2, and if
		#Y1 is bigger than another Y1 and less than the corresponding Y2, then there is
		#a problem
		if(any(bounds$plot.bounds[i,1] > bounds$plot.bounds[j,1] &
			bounds$plot.bounds[i,1] < bounds$plot.bounds[j,2] &
			bounds$plot.bounds[i,3] > bounds$plot.bounds[j,3] &
			bounds$plot.bounds[i,3] < bounds$plot.bounds[j,4]))
		{
			#turn error to TRUE and break the for loop, or else it will get written over
			error[i] <- TRUE
			break;
		}
		else
		{
			error[i] <- FALSE
		}
	}
}

test_that("Plots are sampled and returned in appropriate format",
{
	#cdm should be in matrix format
	expect_is(cdm$picante.cdm, "matrix")
	#plots without any species are cut, so just confirm there are at least some rows
	#species that do not occur are still in cdm, so there should be fifty columns
	expect_true(dim(cdm$picante.cdm)[1] > 1)
	expect_true(dim(cdm$picante.cdm)[2] == 50)
})

test_that("Plots are non-overlapping",
{
	expect_false(any(error))
})
