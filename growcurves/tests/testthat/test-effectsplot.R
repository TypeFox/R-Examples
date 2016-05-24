## library(growcurves, quietly = TRUE)

context("effectsplot returns valid plot object and associated data.frame")

##
## Load simulation dataset without nuisance covariates
## (Two treatment levels, {0,1}, and no nuisance covariates)
##
data(datsimmult)

##
## functions to produce both dpgrowmm and dpgrowmult objects
##
mod.one <- function(x){
		dpgrowmm(y = datsimmult$y, subject = datsimmult$subject, trt = datsimmult$trt, time = datsimmult$time, 
						n.random = datsimmult$n.random, n.fix_degree = 2, Omega = datsimmult$Omega[[1]], group = datsimmult$group[[1]],
						subj.aff = datsimmult$subj.aff[[2]], W.subj.aff = datsimmult$W.subj.aff[[2]], n.iter = 8, n.burn = 2,
					      	n.thin = 1, shape.dp = 4, strength.mm = 0.01, plot.out = TRUE, option = x)  
	 }
	
mod.two <- function(x, niter, nburn, nthin){
		 dpgrowmult(y = datsimmult$y, subject = datsimmult$subject, trt = datsimmult$trt, time = datsimmult$time, 
						n.random = datsimmult$n.random, n.fix_degree = 2, Omega = datsimmult$Omega, group = datsimmult$group,
						subj.aff = datsimmult$subj.aff, W.subj.aff = datsimmult$W.subj.aff, n.iter = 8, n.burn = 2,
					        n.thin = 1, shape.dp = 4, plot.out = TRUE, option = x, ulabs = paste("cbt",1:4,sep="_")) 
	 }



test_that("effectsplot returns plot object when plotting specific MM terms", {

	ONE	<- mod.one("mmcar")
	TWO	<- mod.two(c("mmi","mmcar","mmi","mmcar"))

	run.objects 	= list(ONE = ONE, TWO = TWO)
	mm.terms 	= c(1,2) ## the first object employs the second group (18 clients x 40 sessions) of 4, total, from datsimmult.  There is only one group, so all sessions are selected.
        prior.labs	= c("mmcar","mmcar")
        axis.labs	= c("sessions attended","session effect")

	ep <- effectsplot(objects = run.objects, mm.terms = mm.terms, prior.labs = prior.labs, axis.labs = axis.labs)

	expect_that(length(ep), equals(2)) ## normally would be 4, but here each object has different number of terms, so will only make one plot
	expect_that(ep$p.term, is_a("ggplot"))
	expect_that(names(ep$dat.term)[1], matches("session"))
	expect_that(nrow(ep$dat.term),is_equivalent_to(length(run.objects)*nrow(summary(ONE)$summary.results$u.summary))) ## plot data.frame properly dimensioned

})


test_that("effectsplot returns plot object when plotting ALL terms", {

	TWO	<- mod.two(c("mmi","mmcar","mmi","mmcar"))
	THREE	<- mod.two(c("mmi","mmi","mmi","mmi"))

	run.objects 	= list(TWO = TWO, THREE = THREE)
	mm.terms 	= c("cbt_2","cbt_2")
        prior.labs	= c("mixed","mmi")
        axis.labs	= c("sessions attended","session effect")

	ep <- effectsplot(objects = run.objects, mm.terms = mm.terms, prior.labs = prior.labs, axis.labs = axis.labs)

	expect_that(length(ep), equals(4)) ## both plotting a term and the objects have the same number of  terms
	expect_that(ep$p.all, is_a("ggplot"))
	expect_that(names(ep$dat.all)[1], matches("prior"))

})