## library(growcurves, quietly = TRUE)

context("dpgrowmm returns correct objects")

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
						subj.aff = datsim$subj.aff, W.subj.aff = datsim$W.subj.aff, n.iter = niter, n.burn = nburn,
					      n.thin = nthin, shape.dp = 4, strength.mm = 0.001, plot.out = TRUE, option = x) 
	 }

test_that("mmcar option of dpgrowmm returns expected objects", {
	
	niter		<- 8
	nburn		<- 2
	nthin		<- 2
	MMCAR		<- mod("mmcar",niter,nburn,nthin)
	srm 		<- summary(MMCAR)$summary.results
	parms		<- samples(MMCAR)
	pr		<- MMCAR$plot.results
	num.subj	<- length(unique(datsim$subject))
	nrandom		<- ncol(srm$Z)
	nsessions	<- nrow(datsim$Omega)

	## evaluating class
	expect_that(MMCAR, is_a("dpgrowmm"))
	## evaluating summary output
	expect_that(length(names(srm)), equals(25))
	expect_that(nrow(srm$u.summary), is_equivalent_to(nsessions))
	expect_that(srm$bmat.summary, is_a("list"))
	expect_that(names(srm)[16], matches("lpml"))
	expect_that(ncol(srm$X),equals(5))
	expect_match(colnames(srm$X),"time")
	## evaluating MCMC sample results
	expect_that(nrow(parms$M),is_equivalent_to((niter-nburn)/nthin))
	expect_that(ncol(parms$B),is_equivalent_to(num.subj*nrandom))
	expect_that(length(residuals(MMCAR)),equals(length(datsim$y)))
	## checking plot output
	expect_that(length(names(pr)),equals(11))
	expect_that(names(pr)[11],matches("p.gcsel"))
})


test_that("mmigrp option of dpgrowmm returns expect objects", {
	
	niter		<- 8
	nburn		<- 2
	nthin		<- 2
	MMIGRP		<- mod("mmigrp",niter,nburn,nthin)
	srm 		<- summary(MMIGRP)$summary.results
	parms		<- samples(MMIGRP)
	pr		<- MMIGRP$plot.results
	num.subj	<- length(unique(datsim$subject))
	nrandom		<- ncol(srm$Z)
	nsessions	<- nrow(datsim$Omega)

	## evaluating class
	expect_that(MMIGRP,is_a("dpgrowmm"))
	## evaluating summary output
	expect_that(length(names(srm)), equals(25))
	expect_that(nrow(srm$u.summary), is_equivalent_to(nsessions))
	expect_that(srm$bmat.summary, is_a("list"))
	expect_that(names(srm)[16], matches("lpml"))
	expect_that(ncol(srm$X),equals(5))
	expect_match(colnames(srm$X),"time")
	## evaluating MCMC sample results
	expect_that(nrow(parms$M),is_equivalent_to((niter-nburn)/nthin))
	expect_that(ncol(parms$B),is_equivalent_to(num.subj*nrandom))
	expect_that(length(residuals(MMIGRP)),equals(length(datsim$y)))
	## checking plot output
	expect_that(length(names(pr)),equals(11))
	expect_that(names(pr)[11],matches("p.gcsel"))
	expect_that(pr$p.gcsel,is_a("ggplot"))
})


test_that("mmi option of dpgrowmm returns expect objects", {
	
	niter		<- 8
	nburn		<- 2
	nthin		<- 2
	MMI		<- mod("mmi",niter,nburn,nthin)
	srm 		<- summary(MMI)$summary.results
	parms		<- samples(MMI)
	pr		<- MMI$plot.results
	num.subj	<- length(unique(datsim$subject))
	nrandom		<- ncol(srm$Z)
	nsessions	<- nrow(datsim$Omega)

	## evaluating class
	expect_that(MMI,is_a("dpgrowmm"))
	## evaluating summary output
	expect_that(length(names(srm)), equals(25))
	expect_that(nrow(srm$u.summary), is_equivalent_to(nsessions))
	expect_that(srm$bmat.summary, is_a("list"))
	expect_that(names(srm)[16], matches("lpml"))
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

