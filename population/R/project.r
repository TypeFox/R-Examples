#########################################################################
# MONTECARLO SINGLE CORE
#########################################################################

project <- function(
	years,
	runs,
	initial_population,
	survival,
	litter,
	seed
) {

	if(missing(seed)) seed <- 1
	set.seed(seed)

	classes <- length(initial_population)

	output <- .Call("C_montecarlo",
									as.integer(seed),
									as.integer(years),
									as.integer(runs),
									as.integer(initial_population),
									as.double(survival),
									as.double(litter)
	)

	results <- list()
	results$runs <- aperm(`dim<-`(t(output$runs), list(classes, years, runs)), c(2, 1, 3))

	x <- gc(verbose = FALSE)

	return(results)

}

