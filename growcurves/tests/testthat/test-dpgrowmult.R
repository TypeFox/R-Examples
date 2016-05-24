## library(growcurves, quietly = TRUE)

context("dpgrowmult returns correct objects")

##
## Load simulation dataset without nuisance covariates
## (Two treatment levels, {0,1}, and no nuisance covariates)
##
data(datsimmult)

##
## function to run either mmcar, mmigrp or mmi options under dpgrow function
##
mod <- function(x, niter, nburn, nthin){
		 dpgrowmult(y = datsimmult$y, subject = datsimmult$subject, trt = datsimmult$trt, time = datsimmult$time, 
						n.random = datsimmult$n.random, n.fix_degree = 2, Omega = datsimmult$Omega, group = datsimmult$group,
						subj.aff = datsimmult$subj.aff, W.subj.aff = datsimmult$W.subj.aff, n.iter = niter, n.burn = nburn,
					        n.thin = nthin, strength.mm = 0.1, shape.dp = 4, plot.out = TRUE, option = x, ulabs = paste("cbt",1:4,sep="_")) 
	 }

test_that("dpgrowmult with multiple MM terms returns expected objects", {
	
	option		<- c("mmi","mmcar","mmi","mmcar")
	niter		<- 5
	nburn		<- 2
	nthin		<- 1
	MMMULT		<- mod(option,niter,nburn,nthin)
	srm 		<- summary(MMMULT)$summary.results
	parms		<- samples(MMMULT)
	pr		<- MMMULT$plot.results
	num.subj	<- length(unique(datsimmult$subject))
	nrandom		<- ncol(srm$Z)
	numt		<- srm$numt
	Nterms		<- length(numt)
	option		<- srm$model

	## evaluating class
	expect_that(MMMULT,is_a("dpgrowmult"))
	## evaluating summary output
	## expect_that(length(names(srm)), equals(22))
	expect_that(srm$u.summary, is_a("list"))
	expect_that(length(option), equals(Nterms))
	expect_that(nrow(srm$u.summary[[2]]), is_equivalent_to(numt[2]))
	expect_that(nrow(srm$tauu.summary), is_equivalent_to(Nterms))
	expect_that(srm$bmat.summary, is_a("list"))
	## expect_that(names(srm)[16], matches("lpml"))
	expect_that(ncol(srm$X),equals(5))
	expect_match(colnames(srm$X),"time")
	## evaluating MCMC sample results
	expect_that(length(parms$Gamma), equals(length(option)))
	expect_that(ncol(parms$Gamma[[2]]), is_equivalent_to(numt[2]))
	expect_that(names(parms$Gamma), is_equivalent_to(as.character(srm$ulabs)))
	expect_that(nrow(parms$M),is_equivalent_to((niter-nburn)/nthin))
	expect_that(ncol(parms$B),is_equivalent_to(num.subj*nrandom))
	expect_that(length(residuals(MMMULT)),equals(length(datsimmult$y)))
	## checking plot output
	## expect_that(length(names(pr)),equals(11))
	## expect_that(names(pr)[11],matches("p.gcsel"))
})


