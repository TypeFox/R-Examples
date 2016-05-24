## library(growcurves, quietly = TRUE)

context("ddpgrow returns correct objects")

##
## Load simulation dataset without nuisance covariates
## (Two treatment levels, {0,1}, and no nuisance covariates)
##
data(datsimmult)

##
## function to run dprgrow function under all possible prior choices - c("mvn","car","ind")
##
mod <- function(niter, nburn, nthin, typetreat){
		ddpgrow(y = datsimmult$y, subject = datsimmult$subject, trt = datsimmult$trt, time = datsimmult$time, 
						n.random = datsimmult$n.random, n.fix_degree = 2, dosemat = datsimmult$dosemat_test, Omega = datsimmult$Omega_test,
						numdose = datsimmult$numdose_test, labt = datsimmult$labt_test, typetreat = typetreat, n.iter = niter, n.burn = nburn,
					      	n.thin = nthin, shape.dp = 1.0, rate.dp = 1.0, M.init = length(unique(datsimmult$subject)), plot.out = TRUE)  
	 }

test_that("ddp function returns correct objects under car prior construction", {
	
	niter		<- 30
	nburn		<- 10
	nthin		<- 2
	typetreat	<- c("car")
	DDP		<- mod(niter,nburn,nthin,typetreat)
	srm 		<- summary(DDP)$summary.results
	parms		<- samples(DDP)
	pr		<- DDP$plot.results
	num.subj	<- length(unique(datsimmult$subject))
	nrandom		<- ncol(srm$Z)

	## evaluating class
	expect_that(DDP,is_a("ddpgrow"))
	## evaluating summary output
	expect_that(srm$theta.summary, is_a("list"))
	expect_that(names(srm)[17], matches("lpml"))
	expect_that(ncol(srm$X),equals(5))
	expect_that(colnames(srm$X),matches("time"))
	## evaluating MCMC sample results
	expect_that(nrow(parms$M),is_equivalent_to((niter-nburn)/nthin))
	expect_that(ncol(parms$Theta),is_equivalent_to(num.subj*nrandom))
	expect_that(length(residuals(DDP)),equals(length(datsimmult$y)))
	expect_that(nrow(srm$taucar.summary),equals(1)) ## only 1 car term
	## checking plot output
	expect_that("p.tcar" %in% names(pr), is_true())
	## checking supplemental plot
	psupp	<- ddpEffectsplot(DDP)
	expect_that(length(psupp), equals(4))
})

test_that("ddp function returns correct objects under mvn prior construction", {
	
	niter		<- 30
	nburn		<- 10
	nthin		<- 2
	typetreat	<- c("mvn")
	DDP		<- mod(niter,nburn,nthin,typetreat)
	srm 		<- summary(DDP)$summary.results
	parms		<- samples(DDP)
	pr		<- DDP$plot.results
	num.subj	<- length(unique(datsimmult$subject))
	nrandom	<- ncol(srm$Z)

	## evaluating class
	expect_that(DDP,is_a("ddpgrow"))
	## evaluating summary output
	expect_that(srm$theta.summary, is_a("list"))
	expect_that(names(srm)[17], matches("lpml"))
	expect_that(ncol(srm$X),equals(5))
	expect_that(colnames(srm$X),matches("time"))
	## evaluating MCMC sample results
	expect_that(nrow(parms$M),is_equivalent_to((niter-nburn)/nthin))
	expect_that(ncol(parms$Theta),is_equivalent_to(num.subj*nrandom))
	expect_that(length(residuals(DDP)),equals(length(datsimmult$y)))
	expect_that(length(srm$pmvn.summary),equals(1)) ## only 1 mvn term
	## checking plot output
	expect_that("p.mvn" %in% names(pr), is_true())
	## checking supplemental plot
	psupp	<- ddpEffectsplot(DDP, trts.plot = "num_of_sessions")
	expect_that(length(psupp), equals(4))
})


test_that("ddp function returns correct objects under ind prior construction", {
	
	niter		<- 30
	nburn		<- 10
	nthin		<- 2
	typetreat	<- c("ind")
	DDP		<- mod(niter,nburn,nthin,typetreat)
	srm 		<- summary(DDP)$summary.results
	parms		<- samples(DDP)
	pr		<- DDP$plot.results
	num.subj	<- length(unique(datsimmult$subject))
	nrandom	<- ncol(srm$Z)

	## evaluating class
	expect_that(DDP,is_a("ddpgrow"))
	## evaluating summary output
	expect_that(srm$theta.summary, is_a("list"))
	expect_that(names(srm)[17], matches("lpml"))
	expect_that(ncol(srm$X),equals(5))
	expect_that(colnames(srm$X),matches("time"))
	## evaluating MCMC sample results
	expect_that(nrow(parms$M),is_equivalent_to((niter-nburn)/nthin))
	expect_that(ncol(parms$Theta),is_equivalent_to(num.subj*nrandom))
	expect_that(length(residuals(DDP)),equals(length(datsimmult$y)))
	expect_that(length(srm$tauind.summary),equals(1)) ## only 1 ind term
	## checking plot output
	expect_that("p.iband" %in% names(pr), is_true())
	## checking supplemental plot
	psupp	<- ddpEffectsplot(DDP, trts.plot = "num_of_sessions")
	expect_that(length(psupp), equals(4))
})
