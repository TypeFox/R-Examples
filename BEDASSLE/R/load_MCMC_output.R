load_MCMC_output <- 
function(MCMC.output){
    tmpenv <- environment()
	tmp <- load(MCMC.output,envir=tmpenv)
	mcmc.output <- lapply(tmp,get,envir=tmpenv)
	names(mcmc.output) <- tmp
	stopifnot(length(intersect(names(mcmc.output),c("a0","aD","aE","a2",
									"beta","last.params",
									"LnL_thetas","LnL_counts",
									"Prob","samplefreq","ngen",
									"a0_moves","aD_moves",
									"aE_moves","a2_moves",
									"thetas_moves","mu_moves","beta_moves",
									"aD_accept","aE_accept","a2_accept",
									"thetas_accept","mu_accept")
						)) > 20)
	return(mcmc.output)
}