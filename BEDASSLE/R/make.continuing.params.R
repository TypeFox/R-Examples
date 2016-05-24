make.continuing.params <-
function(MCMC.output,file.name){
	MCMC.output.list <- load_MCMC_output(MCMC.output)
	with(MCMC.output.list, {
		if(!exists("phi_mat")){
			continuing.params <- list(last.params$a0,
										last.params$aD,
										last.params$aE,
										last.params$a2,
										last.params$beta,
										last.params$mu,
										last.params$thetas)
			names(continuing.params) <- c("a0","aD","aE","a2","beta","mu","thetas")										
		}
		if(exists("phi_mat")){
			continuing.params <- list(last.params$a0,
										last.params$aD,
										last.params$aE,
										last.params$a2,
										last.params$beta,
										last.params$phi,
										last.params$mu,
										last.params$thetas)
			names(continuing.params) <- c("a0","aD","aE","a2","beta","phi","mu","thetas")
		}
	save(continuing.params,file=file.name)
	})
}
