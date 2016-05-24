## library(growcurves, quietly = TRUE)

context("trtplot returns valid plot object and associated data.frame")

test_that("trtplot returns plot object", {
	##
	## Load simulation dataset without nuisance covariates
	## (Two treatment levels, {0,1}, and no nuisance covariates)
	##
	data(datsim)

	##
	## function to run either dp or lgm options under dpgrow function
	##
	mod <- function(x){
		dpgrow(y = datsim$y, subject = datsim$subject, trt = datsim$trt, time = datsim$time, 
						n.random = datsim$n.random, n.fix_degree = 2, n.iter = 12, n.burn = 2,
					      n.thin = 2, shape.dp = 1, plot.out = TRUE, option = x) 
	 }

	
	DP	<- mod("dp")
	LGM	<- mod("lgm")

	run.objects = list(DP = DP, LGM = LGM)
	run.models = c("dp","lgm")
	trt.labs = c(0,1)
	time.points = c(3,6)

	tp <- trtplot(run.objects, run.models, trt.labs, time.points, n.thin = 1) 

	expect_that(length(tp), equals(2))
	expect_that(tp$p.trt, is_a("ggplot"))
	expect_that(names(tp$dat.trt)[1], matches("Mu_diff"))
	expect_that(nrow(tp$dat.trt),is_equivalent_to(length(run.models)*length(time.points)*nrow(samples(DP)$Beta))) ## plot data.frame properly dimensioned

})