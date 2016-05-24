#########################################################################
# MONTECARLO MULTI-CORE
#########################################################################

project_parallel <- function(
	years,
	runs,
	initial_population,
	survival,
	litter,
	cores
) {

	if (missing(cores)) {
		cores <- get_cores(runs)
	}

	if (Sys.info()[['sysname']] == 'Windows') {
		cores <- 1
	}

	runs <- floor(runs/cores)

	simulation1 <- mclapply(1:cores,
													project,
													years = years,
													runs = runs,
													initial_population = initial_population,
													survival = survival,
													litter = litter,
													mc.cores = cores
	)

	simulation2 <- list()

	x <- NULL
	for (i in 1:length(simulation1)) {
		x <- abind(x, simulation1[[i]]$runs, along=3)
	}
	simulation2$runs <- x

	return(simulation2)

}
