library(metricTester)
context("Ensure post-simulation analysis functions work if concat.by equals both")

#create a list of results like what you'd expect if concat.by = "both"

ses.null1.byrichness <- matrix(nrow=10, ncol=4)

ses.null1.byrichness <- as.data.frame(ses.null1.byrichness)

ses.null1.byrichness[,1] <- 10:19
set.seed(0)
ses.null1.byrichness[,2] <- rnorm(n=10, mean=0, sd=0.5)
set.seed(0)
ses.null1.byrichness[,3] <- rnorm(n=10, mean=1115, sd=0.5)
set.seed(0)
ses.null1.byrichness[,4] <- rnorm(n=10, mean=-1115, sd=0.5)

names(ses.null1.byrichness) <- c("richness","metric1","metric2","metric3")

ses.null1.byplot <- ses.null1.byrichness
ses.null1.byplot[,1] <- paste("plot", 10:19, sep="")
names(ses.null1.byplot) <- c("plot","metric1","metric2","metric3")

plot.null1.byrichness <- ses.null1.byrichness

plot.null1.byrichness[,2] <- sample(c(0,1,2), 10, replace=TRUE)
plot.null1.byrichness[,3] <- sample(c(0,1,2), 10, replace=TRUE)
plot.null1.byrichness[,4] <- sample(c(0,1,2), 10, replace=TRUE)

plot.null1.byplot <- plot.null1.byrichness
plot.null1.byplot[,1] <- paste("plot", 10:19, sep="")
names(plot.null1.byplot) <- c("plot","metric1","metric2","metric3")

ses.null1 <- list("richness"=ses.null1.byrichness, "plot"=ses.null1.byplot)

ses1 <- list("null1"=ses.null1, "null2"= ses.null1)

plot.null1 <- list("richness"=plot.null1.byrichness,
	"plot"=plot.null1.byplot)

plot1 <- list("null1"=plot.null1, "null2"=plot.null1)

random1 <- list("ses"=ses1, "plot"=plot1)

fake.results1 <- list("random"=random1, "filtering"=random1, "competition"=random1)

results <- list("iteration1"=fake.results1, "iteration2"=fake.results1)

#now summarize these results
resultsSumm <- reduceResults(results, "both")

#generate some temporary results
failedTemp <- failed(results[[1]], "both")
sesSingleTemp <- sesSingle(results[[1]], "both")
sesIndivTemp <- sesIndiv(results, "both")
sesOverallTemp <- sesOverall(resultsSumm$ses, test="wilcotest", concat.by="both")

#confirm the failed function returns a data frame with no rows and four columns
#when provided with results concatenated by both
test_that("data frame of correct dimensions is returned",
{
	expect_is(failedTemp, "data.frame")
	expect_true(dim(failedTemp)[1] == 0)
	expect_true(dim(failedTemp)[2] == 4)
})

#confirm that sesSingle works on a single iteration results concatenated by both
test_that("data frame of correct dimensions is returned",
{
	expect_is(sesSingleTemp, "data.frame")
	expect_true(dim(sesSingleTemp)[1] == 36)
	expect_true(dim(sesSingleTemp)[2] == 6)
})

test_that("sesSingle is correctly naming the summarized results",
{
	expect_true(all(sesSingleTemp[,"simulation"] == c(rep("random",12),
		rep("filtering", 12), rep("competition", 12))))
	expect_true(all(sesSingleTemp[,"null.model"] == rep(c(rep("null1",6),
		rep("null2", 6)),6)))
	expect_true(all(sesSingleTemp[,"concat.by"] == rep(c(rep("richness", 3),
		rep("plot",3)),6)))
	expect_true(all(sesSingleTemp[,"metric"] == rep(c("metric1", "metric2","metric3"),
		12)))
})

test_that("p values make sense for the different spatial simulations",
{
	#if random and SES scores are centered around zero, p should be large
	expect_true(sesSingleTemp[1,"p.value"] > 0.05)
	#if random and SES scores not centered around zero, p should be small
	expect_true(sesSingleTemp[2,"p.value"] <= 0.05)
	expect_true(sesSingleTemp[3,"p.value"] <= 0.05)
	#if filtering and SES scores are centered around zero or around a positive value
	#p should be large. should be small if centered around a negative value
	expect_true(sesSingleTemp[13,"p.value"] > 0.05)
	expect_true(sesSingleTemp[14,"p.value"] > 0.05)
	expect_true(sesSingleTemp[15,"p.value"] <= 0.05)	
	#if competition and SES scores are centered around zero or around a negative value
	#p should be large. should be small if centered around a positive value
	expect_true(sesSingleTemp[25,"p.value"] > 0.05)
	expect_true(sesSingleTemp[26,"p.value"] <= 0.05)
	expect_true(sesSingleTemp[27,"p.value"] > 0.05)	
})

#confirm that sesIndiv works on multiple iteration results concatenated by both
test_that("data frame of correct dimensions is returned",
{
	expect_is(sesIndivTemp, "data.frame")
	expect_true(dim(sesIndivTemp)[1] == 36)
	expect_true(dim(sesIndivTemp)[2] == 7)
})

test_that("sesIndiv is correctly naming the summarized results",
{
	expect_true(all(sesIndivTemp[,"simulation"] == c(rep("random",12),
		rep("filtering", 12), rep("competition", 12))))
	expect_true(all(sesIndivTemp[,"null.model"] == rep(c(rep("null1",6),
		rep("null2", 6)),6)))
	expect_true(all(sesIndivTemp[,"concat.by"] == rep(c(rep("richness", 3),
		rep("plot",3)),6)))
	expect_true(all(sesIndivTemp[,"metric"] == rep(c("metric1", "metric2","metric3"),
		12)))
})

test_that("sesIndiv output makes sense for the different spatial simulations",
{
	#expecting two runs of each, so
	expect_true(sum(sesIndivTemp[,"total.runs"]) == 72)
	#both runs should have not thrown a type I error for metric 1 & should have for
	#metrics 2 & 3 for random spatial sim.
	expect_true(sesIndivTemp[1,"typeI"] == 0)
	expect_true(sesIndivTemp[2,"typeI"] == 2)
	expect_true(sesIndivTemp[3,"typeI"] == 2)
	#type II for random spatial sims should be NA
	expect_true(is.na(sesIndivTemp[1,"typeII"]))
	#both runs should throw typeII errors for metrics 1 & 2 and not for 3 for filtering
	expect_true(sesIndivTemp[13,"typeII"] == 2)
	expect_true(sesIndivTemp[14,"typeII"] == 2)
	expect_true(sesIndivTemp[15,"typeII"] == 0)
	#type I for filtering spatial sims should be NA
	expect_true(is.na(sesIndivTemp[13,"typeI"]))
	#both runs should throw typeII errors for metrics 1 & 3 and not for 2 for competition
	expect_true(sesIndivTemp[25,"typeII"] == 2)
	expect_true(sesIndivTemp[26,"typeII"] == 0)
	expect_true(sesIndivTemp[27,"typeII"] == 2)
	#type I for filtering spatial sims should be NA
	expect_true(is.na(sesIndivTemp[25,"typeI"]))
})

#confirm that sesOverall works on multiple iteration results concatenated by both
test_that("data frame of correct dimensions is returned",
{
	expect_is(sesOverallTemp, "data.frame")
	expect_true(dim(sesOverallTemp)[1] == 36)
	expect_true(dim(sesOverallTemp)[2] == 6)
})

test_that("sesIndiv output makes sense for the different spatial simulations",
{
	#both runs should have not thrown a type I error for metric 1 & should have for
	#metrics 2 & 3 for random spatial sim.
	expect_true(sesOverallTemp[1,"p.value"] > 0.01)
	expect_true(sesOverallTemp[2,"p.value"] <= 0.05)
	expect_true(sesOverallTemp[3,"p.value"] <= 0.05)
	#both runs should throw typeII errors for metrics 1 & 2 and not for 3 for filtering
	expect_true(sesOverallTemp[13,"p.value"] > 0.01)
	expect_true(sesOverallTemp[14,"p.value"] > 0.05)
	expect_true(sesOverallTemp[15,"p.value"] <= 0.05)
	#both runs should throw typeII errors for metrics 1 & 3 and not for 2 for competition
	expect_true(sesOverallTemp[25,"p.value"] > 0.01)
	expect_true(sesOverallTemp[26,"p.value"] <= 0.05)
	expect_true(sesOverallTemp[27,"p.value"] > 0.05)
})
