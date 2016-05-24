library(metricTester)
context("Ensure post-simulation analysis functions work if concat.by equals richness")

#create a list of results like what you'd expect if concat.by = "richness"

ses.null1 <- matrix(nrow=10, ncol=4)

ses.null1 <- as.data.frame(ses.null1)

ses.null1[,1] <- 10:19
set.seed(0)
ses.null1[,2] <- rnorm(n=10, mean=0, sd=0.5)
set.seed(0)
ses.null1[,3] <- rnorm(n=10, mean=1000, sd=0.5)
set.seed(0)
ses.null1[,4] <- rnorm(n=10, mean=-1000, sd=0.5)

names(ses.null1) <- c("richness","metric1","metric2","metric3")

plot.null1 <- ses.null1

plot.null1[,2] <- sample(c(0,1,2), 10, replace=TRUE)
plot.null1[,3] <- sample(c(0,1,2), 10, replace=TRUE)
plot.null1[,4] <- sample(c(0,1,2), 10, replace=TRUE)

ses.null2 <- ses.null1
plot.null2 <- plot.null1

ses1 <- list("null1"=ses.null1, "null2"=ses.null2)
plot1 <- list("null1"=plot.null1, "null2"=plot.null2)

random1 <- list("ses"=ses1, "plot"=plot1)

fake.results1 <- list("random"=random1, "filtering"=random1, "competition"=random1)

results <- list("iteration1"=fake.results1, "iteration2"=fake.results1)

#now summarize these results
resultsSumm <- reduceResults(results, "richness")

#generate some temporary results
failedTemp <- failed(results[[1]], "richness")
sesSingleTemp <- sesSingle(results[[1]], "richness")
sesIndivTemp <- sesIndiv(results, "richness")
sesOverallTemp <- sesOverall(resultsSumm$ses, test="wilcotest", concat.by="richness")

#confirm the failed function returns a data frame with no rows and four columns
#when provided with results concatenated by richness
test_that("data frame of correct dimensions is returned",
{
	expect_is(failedTemp, "data.frame")
	expect_true(dim(failedTemp)[1] == 0)
	expect_true(dim(failedTemp)[2] == 4)
})

#confirm that sesSingle works on a single iteration results concatenated by richness
test_that("data frame of correct dimensions is returned",
{
	expect_is(sesSingleTemp, "data.frame")
	expect_true(dim(sesSingleTemp)[1] == 18)
	expect_true(dim(sesSingleTemp)[2] == 6)
})

test_that("sesSingle is correctly naming the summarized results",
{
	expect_true(all(sesSingleTemp[,"simulation"] == c(rep("random",6), rep("filtering", 6),
		rep("competition", 6))))
	expect_true(all(sesSingleTemp[,"null.model"] == rep(c(rep("null1",3),
		rep("null2", 3)),3)))
	expect_true(all(sesSingleTemp[,"concat.by"] == rep("richness", 18)))
	expect_true(all(sesSingleTemp[,"metric"] == rep(c("metric1", "metric2","metric3"), 6)))
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
	expect_true(sesSingleTemp[7,"p.value"] > 0.05)
	expect_true(sesSingleTemp[8,"p.value"] > 0.05)
	expect_true(sesSingleTemp[9,"p.value"] <= 0.05)	
	#if competition and SES scores are centered around zero or around a negative value
	#p should be large. should be small if centered around a positive value
	expect_true(sesSingleTemp[13,"p.value"] > 0.05)
	expect_true(sesSingleTemp[14,"p.value"] <= 0.05)
	expect_true(sesSingleTemp[15,"p.value"] > 0.05)	
})

#confirm that sesIndiv works on multiple iteration results concatenated by richness
test_that("data frame of correct dimensions is returned",
{
	expect_is(sesIndivTemp, "data.frame")
	expect_true(dim(sesIndivTemp)[1] == 18)
	expect_true(dim(sesIndivTemp)[2] == 7)
})

test_that("sesIndiv is correctly naming the summarized results",
{
	expect_true(all(sesIndivTemp[,"simulation"] == c(rep("random",6), rep("filtering", 6),
		rep("competition", 6))))
	expect_true(all(sesIndivTemp[,"null.model"] == rep(c(rep("null1",3),
		rep("null2", 3)),3)))
	expect_true(all(sesIndivTemp[,"concat.by"] == rep("richness", 18)))
	expect_true(all(sesIndivTemp[,"metric"] == rep(c("metric1", "metric2","metric3"), 6)))
})

test_that("sesIndiv output makes sense for the different spatial simulations",
{
	#expecting two runs of each, so
	expect_true(sum(sesIndivTemp[,"total.runs"]) == 36)
	#both runs should have not thrown a type I error for metric 1 & should have for
	#metrics 2 & 3 for random spatial sim.
	expect_true(sesIndivTemp[1,"typeI"] == 0)
	expect_true(sesIndivTemp[2,"typeI"] == 2)
	expect_true(sesIndivTemp[3,"typeI"] == 2)
	#type II for random spatial sims should be NA
	expect_true(is.na(sesIndivTemp[1,"typeII"]))
	#both runs should throw typeII errors for metrics 1 & 2 and not for 3 for filtering
	expect_true(sesIndivTemp[7,"typeII"] == 2)
	expect_true(sesIndivTemp[8,"typeII"] == 2)
	expect_true(sesIndivTemp[9,"typeII"] == 0)
	#type I for filtering spatial sims should be NA
	expect_true(is.na(sesIndivTemp[7,"typeI"]))
	#both runs should throw typeII errors for metrics 1 & 3 and not for 2 for competition
	expect_true(sesIndivTemp[13,"typeII"] == 2)
	expect_true(sesIndivTemp[14,"typeII"] == 0)
	expect_true(sesIndivTemp[15,"typeII"] == 2)
	#type I for filtering spatial sims should be NA
	expect_true(is.na(sesIndivTemp[13,"typeI"]))
})

#confirm that sesOverall works on multiple iteration results concatenated by richness
test_that("data frame of correct dimensions is returned",
{
	expect_is(sesOverallTemp, "data.frame")
	expect_true(dim(sesOverallTemp)[1] == 18)
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
	expect_true(sesOverallTemp[7,"p.value"] > 0.01)
	expect_true(sesOverallTemp[8,"p.value"] > 0.05)
	expect_true(sesOverallTemp[9,"p.value"] <= 0.05)
	#both runs should throw typeII errors for metrics 1 & 3 and not for 2 for competition
	expect_true(sesOverallTemp[13,"p.value"] > 0.01)
	expect_true(sesOverallTemp[14,"p.value"] <= 0.05)
	expect_true(sesOverallTemp[15,"p.value"] > 0.05)
})
