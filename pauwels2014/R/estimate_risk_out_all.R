if(getRversion() >= "2.15.1")  utils::globalVariables("observables")
estimate_risk_out_all <- function(thetas, knobj, experiment_fun){
	## Estimates the risk of experiment defined by experiment_fun given
	## sample thetas and parameters in knobj
	## N_simu_weight is the number of sample required for the output
	## simulation


	## Simulate trajectories related to our experiment for each sample point
	simulation_thetas <- lapply(1:nrow(thetas), FUN = function(i){ simulate_experiment(thetas[i,], knobj, experiment_fun) })

	to_delete <- sapply(simulation_thetas, FUN = function(x){sum(is.nan(x)) > 0})
	#if (sum(to_delete) > 0){
	#     thetas <- thetas[!to_delete,]
	#     simulation_thetas <- lapply(1:nrow(thetas), FUN = function(i){ simulate_experiment(thetas[i,], knobj, experiment_fun) })
	#}

	risksRes <- 1:length(observables)

	## For each observable

	for(obs_id in 1:length(observables)){
		obs_exp <- observables[[obs_id]]$obs
		tspan_exp <- observables[[obs_id]]$reso

		simulation_thetas_exp <- lapply(simulation_thetas, FUN = function(x){ x[x[,1] %in% tspan_exp, dimnames(x)[[2]] %in% obs_exp] })
		simulation_thetas_exp <- array(unlist(simulation_thetas_exp), dim = c(dim(simulation_thetas_exp[[1]]), length(simulation_thetas_exp)), dimnames = list(dimnames(simulation_thetas_exp[[1]])[[1]], dimnames(simulation_thetas_exp[[1]])[[2]], NULL) )

		risk_theta_T <- matrix(0, nrow(thetas), 1)
		for( i_T in 1:nrow(thetas)){
			## For each theta, compute the risk if it was the true parameter
			theta_T <- thetas[i_T,]

			data_theta_T <- simulation_thetas_exp[,,i_T]
			data_theta_Ts <- data_theta_T

			weights <- matrix(0, knobj$global_parameters$n_simu_weights, nrow(thetas))
			data_theta_T_arr <- array(0, c(dim(data_theta_Ts), knobj$global_parameters$n_simu_weights)  )

			## Simulate from noise model
			for (i in 1:knobj$global_parameters$n_simu_weights){
				data_theta_T <- add_noise(data_theta_Ts)
				weights[i,] <- log_normalize(.Call("eval_weights_risk", simulation_thetas_exp, data_theta_T))[[1]]
			}

			expectation_weights <- apply(weights, 2, mean)
			risks_theta <- rep(0, nrow(thetas))
			for (i in 1:nrow(thetas)){
				risks_theta[i] <- risk_theta_fun(knobj$transform_params(theta_T), knobj$transform_params(thetas[i,]), knobj$global_parameters$n_params)
			}
			risk_theta_T[i_T] <- sum(expectation_weights * risks_theta)
			#print(i_T)
		}
		risksRes[obs_id] <- mean(risk_theta_T)
	} 

	res <- c()

	res <- data.frame(Measurement = sapply(observables, FUN = function(x){x$name}), Risk = risksRes, Cost = sapply(observables, FUN = function(x){x$cost}))
	res
}
