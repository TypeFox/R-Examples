if(getRversion() >= "2.15.1")  utils::globalVariables(c("experiment_list1", "observables"))
dream6_design <-
function(knobj, sample_function, seed, credits = 5000, file_to_save = NULL, verbose = T){
	# Perform active design simulation
	# Iterates until all credit has been spent
	# sample_function generates a sample of the posterior distribution
	
	next_it <- TRUE
	k <- 1
	
	while(next_it){
		# Sample and estimate risk
		
		if(verbose){
			print(paste("Sample", k))
		}	
		thetas <- sample_function(knobj)	
		knobj$datas[[length(knobj$datas)]]$thetas_est <- thetas
		knobj$datas[[length(knobj$datas)]]$thetas <- thetas[sample(1:nrow(thetas), size = knobj$global_parameters$final_sample_design),]
		risks <- c()
		if(!is.null(file_to_save)){
			saveRDS(knobj, file_to_save)
		}
		if(verbose){
			print(paste("Estimate risk", k))
		}
		for(id_exp in 1:length(experiment_list1)){
			experiment_fun <- experiment_list1[[id_exp]]
			experiment_list1[[id_exp]]
			res <- estimate_risk_dream6(thetas, knobj, experiment_fun)
			res$Cost <- res$Cost + experiment_fun(NULL, NULL)$cost
			res$exp <- names(experiment_list1)[id_exp]
			risks <- rbind(risks,res)
			print(id_exp)
		}
	
	
		# Choose next experiment and generate data
		if(verbose){
			print(paste("Get data", k))
		}
		knobj$datas[[length(knobj$datas)]]$risks <- risks
		if(!is.null(file_to_save)){
			saveRDS(knobj, file_to_save)
		}
		risks <- risks[risks$Cost <= credits,]
		
		temp_risk <- risks$Risk/risks$Cost
		next_exp <- which.max(temp_risk)[1]
		nnext_it <- paste(risks$exp[next_exp], risks$Measurement[next_exp]) %in% knobj$experiments
		
		## Remove previously seen experiments
		while(nnext_it){
			temp_risk[next_exp] <- 0
			next_exp <- which.max(temp_risk)[1]
			nnext_it <- (paste(risks$exp[next_exp], risks$Measurement[next_exp]) %in% knobj$experiments) & (max(temp_risk)[1] > 0)
		}
		
		if(max(temp_risk)[1] <= 0){
			break
		}
		
		# Update experiment list
		knobj$experiments <- c(knobj$experiments, paste(risks$exp[next_exp], risks$Measurement[next_exp]) )
		exp_fun_next_exp <- experiment_list1[[which(names(experiment_list1) == risks[next_exp, 4])]]
		
		# Simulate data, add noise and take the subset corresponding to the chosen quantity to observe
		data_next_exp <- simulate_experiment(knobj$global_parameters$true_params_T, knobj, exp_fun_next_exp)
		data_next_exp <- add_noise(data_next_exp)		
		to_observe <- observables[[ as.character(risks[next_exp,1]) ]]$obs
		time_res <- observables[[ as.character(risks[next_exp,1]) ]]$reso
		knobj$datas[[length(knobj$datas) + 1]] <- list(manip =  experiment_list1[[which(names(experiment_list1) == risks[next_exp, 4])]], data = data_next_exp[data_next_exp[,1] %in% time_res,to_observe] )
		
		## Pay for the data! and stop if there is not enough to buy new experiments
		credits <- credits - risks[next_exp,]$Cost
		next_it <- (credits >= min(sapply(observables, FUN= function(x){x$cost})))
	
		if(!is.null(file_to_save)){
			saveRDS(knobj, file_to_save)
		}
		k <- k+1
	}
	
	if(verbose){
		print(paste("Sample", k))
	}
	thetas <- sample_function(knobj)	
	knobj$datas[[length(knobj$datas)]]$thetas_est <- thetas
	knobj$datas[[length(knobj$datas)]]$thetas <- thetas[sample(1:nrow(thetas), size = knobj$global_parameters$final_sample_design),]
	if(!is.null(file_to_save)){
		saveRDS(knobj, file_to_save)
	}
	
	knobj
}
