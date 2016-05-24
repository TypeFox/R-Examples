#########################################################################
# MONTECARLO SINGLE CORE
#########################################################################

project <- function(
	years,
	runs,
	surv_pup,
	surv_sub,
	surv_vag,
	surv_adt,
	dispers_weib_shape,
	dispers_weib_scale,
	settl_weib_shape,
	settl_weib_scale,
	pair1breed,
	litter_size,
	pop_initial,
	pop_quota,
	seed
) {

	if(missing(pop_initial)) {
		pop_initial <- list()
		pop_initial$packs <- matrix(rep(c(2,2,5),4), ncol=3, nrow=4, byrow=T)
		pop_initial$vagrants <- 5
	}

	if(missing(pop_quota)) {
		pop_quota <- matrix(0, ncol=12*years+1, nrow=5)
	}

	if(missing(seed)) {
		seed <- 1
	}
	set.seed(seed)

	months <- (years*(12)+1)
	stats <- 15

	output <- .Call("C_montecarlo",
									as.integer(years),
									as.integer(runs),
									as.double(surv_pup),
									as.double(surv_sub),
									as.double(surv_vag),
									as.double(surv_adt),
									as.double(dispers_weib_shape),
									as.double(dispers_weib_scale),
									as.double(settl_weib_shape),
									as.double(settl_weib_scale),
									as.double(pair1breed),
									as.double(litter_size),
									as.integer(as.vector(data.matrix(pop_quota))),
									as.integer(as.vector(data.matrix(pop_initial$packs))),
									as.integer(pop_initial$vagrants)
	)

	results <- list()
	results$runs <- aperm(`dim<-`(t(output$runs), list(stats, months, runs)), c(2, 1, 3))
	results$individuals <- matrix(output$individuals, nrow=length(output$individuals)/5, ncol=5, byrow=T)

	colnames(results$runs) <- c("pop_size",
															"alphas", "alphapairs", "alphagroups", "alphasingle",
															"vagrs", "subs", "pups", "F", "M", "age",
															"pairs", "packs", "psize", "repros")

	colnames(results$individuals) <- c("dispersed", "settled", "firstbred", "died", "run")

	results$parameters <- list(	years = years,
															runs = runs,
															surv_pup = surv_pup,
															surv_sub = surv_sub,
															surv_vag = surv_vag,
															surv_adt = surv_adt,
															dispers_weib_shape = dispers_weib_shape,
															dispers_weib_scale = dispers_weib_scale,
															settl_weib_shape = settl_weib_shape,
															settl_weib_scale = settl_weib_scale,
															pair1breed = pair1breed,
															litter_size = litter_size,
															pop_quota = pop_quota,
															pop_initial = pop_initial)

	x <- gc(verbose = FALSE)

	return(results)

}

