
context("dpgrowmm returns correct objects under multivariate MM effects")

##
## Load simulation dataset without nuisance covariates
## (Two treatment levels, {0,1}, and no nuisance covariates)
##
data(datsim)

##
## function to run either mmcar, mmigrp or mmi options under dpgrow function
##
mod <- function(x, niter, nburn, nthin){
		 dpgrowmm(y = datsim$y, subject = datsim$subject, trt = datsim$trt, time = datsim$time, 
						n.random = datsim$n.random, n.fix_degree = 2, Omega = datsim$Omega, group = datsim$group,
						subj.aff = datsim$subj.aff, W.subj.aff = datsim$W.subj.aff, multi = TRUE, n.iter = niter, n.burn = nburn,
					        n.thin = nthin, shape.dp = 4, strength.mm = 0.1, plot.out = TRUE, option = x) 
	 }


test_that("mmi option of dpgrowmm returns expected objects under multivariate MM effects", {
	
	niter		<- 8
	nburn		<- 2
	nthin		<- 2
	MMI		<- mod("mmi",niter,nburn,nthin)
	srm 		<- summary(MMI)$summary.results
	parms		<- samples(MMI)
	pr		<- MMI$plot.results
	num.subj	<- length(unique(datsim$subject))
	nrandom		<- ncol(srm$Z)
	nsessions	<- ncol(datsim$W.subj.aff)
	Nmv		<- MMI$summary.results$Nmv

	## evaluating class
	expect_that(MMI,is_a("dpgrowmm"))
	## evaluating summary output
	expect_that(length(names(srm)), equals(26))
	expect_that(nrow(srm$u.summary), is_equivalent_to((Nmv*nsessions)))
	expect_that(srm$bmat.summary, is_a("list"))
	expect_that(names(srm)[17], matches("lpml"))
	expect_that(ncol(srm$X),equals(5))
	expect_match(colnames(srm$X),"time")
	## evaluating MCMC sample results
	expect_that(nrow(parms$M),is_equivalent_to((niter-nburn)/nthin))
	expect_that(ncol(parms$B),is_equivalent_to(num.subj*nrandom))
	expect_that(length(residuals(MMI)),equals(length(datsim$y)))
	## checking plot output
	expect_that(length(names(pr)),equals(11))
	expect_that(names(pr)[11],matches("p.gcsel"))
	expect_that(pr$p.gcsel,is_a("ggplot"))
})

test_that("mmcar option of dpgrowmm returns expected objects under multivariate MM effects", {
	
	niter		<- 8
	nburn		<- 2
	nthin		<- 2
	MMCAR		<- mod("mmcar",niter,nburn,nthin)
	srm 		<- summary(MMCAR)$summary.results
	parms		<- samples(MMCAR)
	pr		<- MMCAR$plot.results
	num.subj	<- length(unique(datsim$subject))
	nrandom		<- ncol(srm$Z)
	nsessions	<- ncol(datsim$W.subj.aff)
	Nmv		<- MMCAR$summary.results$Nmv

	## evaluating class
	expect_that(MMCAR,is_a("dpgrowmm"))
	## evaluating summary output
	expect_that(length(names(srm)), equals(26))
	expect_that(nrow(srm$u.summary), is_equivalent_to((Nmv*nsessions)))
	expect_that(srm$bmat.summary, is_a("list"))
	expect_that(names(srm)[17], matches("lpml"))
	expect_that(ncol(srm$X),equals(5))
	expect_match(colnames(srm$X),"time")
	## evaluating MCMC sample results
	expect_that(nrow(parms$M),is_equivalent_to((niter-nburn)/nthin))
	expect_that(ncol(parms$B),is_equivalent_to(num.subj*nrandom))
	expect_that(length(residuals(MMCAR)),equals(length(datsim$y)))
	## checking plot output
	expect_that(length(names(pr)),equals(11))
	expect_that(names(pr)[11],matches("p.gcsel"))
	expect_that(pr$p.gcsel,is_a("ggplot"))
})