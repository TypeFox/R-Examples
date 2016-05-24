## library(growcurves, quietly = TRUE)

context("dpgrow returns correct objects")

##
## Load simulation dataset without nuisance covariates
## (Two treatment levels, {0,1}, and no nuisance covariates)
##
data(datsim)

##
## function to run either dp or lgm options under dpgrow function
##
mod <- function(x, niter, nburn, nthin){
		dpgrow(y = datsim$y, subject = datsim$subject, trt = datsim$trt, time = datsim$time, 
						n.random = datsim$n.random, n.fix_degree = 2, n.iter = niter, n.burn = nburn,
					      n.thin = nthin, shape.dp = 1, plot.out = TRUE, option = x)  
	 }

test_that("dp option of dpgrow returns expect objects", {
	
	niter		<- 8
	nburn		<- 2
	nthin		<- 2
	DP		<- mod("dp",niter,nburn,nthin)
	srm 		<- summary(DP)$summary.results
	parms		<- samples(DP)
	pr		<- DP$plot.results
	num.subj	<- length(unique(datsim$subject))
	nrandom		<- ncol(srm$Z)

	## evaluating class
	expect_that(DP,is_a("dpgrow"))
	## evaluating summary output
	expect_that(length(names(srm)), equals(22))
	expect_that(srm$bmat.summary, is_a("list"))
	expect_that(names(srm)[15], matches("lpml"))
	expect_that(ncol(srm$X),equals(5))
	expect_that(colnames(srm$X),matches("time"))
	## evaluating MCMC sample results
	expect_that(nrow(parms$M),is_equivalent_to((niter-nburn)/nthin))
	expect_that(ncol(parms$B),is_equivalent_to(num.subj*nrandom))
	expect_that(length(residuals(DP)),equals(length(datsim$y)))
	## checking plot output
	expect_that(length(names(pr)),equals(7))
	expect_that(names(pr)[7],matches("p.gcsel"))
})


test_that("lgm option of dpgrow returns expect objects", {
	
	niter		<- 8
	nburn		<- 2
	nthin		<- 2
	LGM		<- mod("lgm",niter,nburn,nthin)
	srm 		<- summary(LGM)$summary.results
	parms		<- samples(LGM)
	pr		<- LGM$plot.results
	num.subj	<- length(unique(datsim$subject))
	nrandom	<- ncol(srm$Z)

	## evaluating class
	expect_that(LGM,is_a("dpgrow"))
	## evaluating summary output
	expect_that(length(names(srm)), equals(22))
	expect_that(srm$bmat.summary, is_a("list"))
	expect_that(names(srm)[15], matches("lpml"))
	expect_that(ncol(srm$X),equals(5))
	expect_that(colnames(srm$X),matches("time"))
	## evaluating MCMC sample results
	expect_that(nrow(parms$Tau.e),is_equivalent_to((niter-nburn)/nthin))
	expect_that(ncol(parms$B),is_equivalent_to(num.subj*nrandom))
	expect_that(length(residuals(LGM)),equals(length(datsim$y)))
	## checking plot output
	expect_that(length(names(pr)),equals(6))
	expect_that(names(pr)[6],matches("p.gcsel"))
})


test_that("dpgrow runs correctly with nuisance fixed effects", {
	
	data(datsimcov)
	niter		<- 8
	nburn		<- 2
	nthin		<- 2
	DP				= dpgrow(y = NULL, subject = datsimcov$subject, trt = datsimcov$trt, time = datsimcov$time,
						n.random = datsimcov$n.random, n.fix_degree = 2, formula = datsimcov$formula, random.only = FALSE, 
						data = datsimcov$data, n.iter = niter, n.burn = nburn, n.thin = nthin, shape.dp = 4, plot.out = TRUE, 
						option = "dp")
	srm 		<- summary(DP)$summary.results
	parms		<- samples(DP)
	pr		<- DP$plot.results
	num.subj	<- length(unique(datsim$subject))
	nrandom	<- ncol(srm$Z)


	## evaluating class
	expect_that(DP,is_a("dpgrow"))
	## evaluating summary output
	expect_that(length(names(srm)), equals(22))
	expect_that(srm$bmat.summary, is_a("list"))
	expect_that(names(srm)[15], matches("lpml"))
	expect_that(ncol(srm$X),equals(7))
	expect_that(colnames(srm$X)[7],matches("income"))
	## evaluating MCMC sample results
	expect_that(nrow(parms$M),is_equivalent_to((niter-nburn)/nthin))
	expect_that(ncol(parms$B),is_equivalent_to(num.subj*nrandom))
	expect_that(length(residuals(DP)),equals(length(datsim$y)))
	## checking plot output
	expect_that(length(names(pr)),equals(7))
	expect_that(names(pr)[7],matches("p.gcsel"))

})


