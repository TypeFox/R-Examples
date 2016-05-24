if(getRversion() >= "2.15.1")  utils::globalVariables("observables")
estimate_risk_dream6 <- function(thetas, knobj, experiment_fun){
	## Estimates the risk of experiment defined by experiment_fun given
	## sample thetas and parameters in knobj
	## N_simu_weight is the number of sample required for the output
	## simulation


	## Simulate trajectories related to our experiment for each sample point
	#simulation_thetas <- lapply(1:nrow(thetas), FUN = function(i){ simulate_experiment(thetas[i,], knobj, experiment_fun) })
	simulation_thetas <- sapply(1:nrow(thetas), FUN = function(i){ simulate_experiment(thetas[i,], knobj, experiment_fun) }, simplify = "array")
	
	

	to_delete <- sapply(simulation_thetas, FUN = function(x){sum(is.nan(x)) > 0})
	#if (sum(to_delete) > 0){
	#     thetas <- thetas[!to_delete,]
	#     simulation_thetas <- lapply(1:nrow(thetas), FUN = function(i){ simulate_experiment(thetas[i,], knobj, experiment_fun) })
	#}
	
	temp_score <- apply(simulation_thetas,c(1,2), function(x){temp <- sqrt(var(x)); if(temp == 0){temp <- 1}; (max(x) - min(x))/temp})

	risksRes <- 1:length(observables)

	## For each observable

	for(obs_id in 1:length(observables)){
		obs_exp <- observables[[obs_id]]$obs
		tspan_exp <- observables[[obs_id]]$reso

		#simulation_thetas_exp <- lapply(simulation_thetas, FUN = function(x){ x[x[,1] %in% tspan_exp, dimnames(x)[[2]] %in% obs_exp] })
		scores <- temp_score[simulation_thetas[,1,1] %in% tspan_exp, dimnames(temp_score)[[2]] %in% obs_exp] 
		
		risksRes[obs_id] <- sum(apply(scores,2,max)) / (ncol(scores) - 1)
	} 

	res <- c()

	res <- data.frame(Measurement = sapply(observables, FUN = function(x){x$name}), Risk = risksRes, Cost = sapply(observables, FUN = function(x){x$cost}))
	res
}
