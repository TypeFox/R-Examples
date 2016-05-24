## library(growcurves, quietly = TRUE)

context("trtplot returns valid plot object and associated data.frame")

test_that("growplot returns plot object", {

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
					        n.thin = 1, shape.dp = 1, plot.out = TRUE, option = x) 
	 }
	
	DP	<- mod("dp")
   	LGM	<- mod("lgm")
	subjects.plot = sample(datsim$subject,10,replace = FALSE)

	gp <- growplot(object = DP, compare.objects = list(LGM = LGM), subjects.plot = subjects.plot, main.label = "DP")

	## determine how many repeated subject measures associated with subjects in subjects.plot
	tmp 			<- data.frame(datsim$y,datsim$subject)
	names(tmp) 		<- c("y","subject") 
	tmp			<- subset(tmp, subject %in% subjects.plot)
	num.y			<- length(tmp$y)
	rm(tmp)

	expect_that(length(gp), equals(5))
	expect_that(gp$p.gctrt, is_a("ggplot"))
	expect_that(names(gp$dat.gc)[1], matches("fit"))
	expect_that(nrow(gp$dat.data),is_equivalent_to(num.y)) ## plot data.frame properly dimensioned

})