if(getRversion() >= "2.15.1")  utils::globalVariables(c("experiment_list1", "observables"))
random_design <-
function(knobj, sample_function, exps, seed, credits = 5000, file_to_save = NULL, verbose = T){
	## Perform random design simulation
	## Iterates until all credit has been spent
	## sample_function generates a sample of the posterior distribution
	## exps shows combination of measurements and perturbations with their cost
	
	next_it <- TRUE
	k <- 1
	
	
	while(next_it){
		## Sample and estimate risk
		
		if(verbose){
			print(paste("Sample", k))
		}
		thetas <- sample_function(knobj)	
		knobj$datas[[length(knobj$datas)]]$thetas_est <- thetas
		knobj$datas[[length(knobj$datas)]]$thetas <- thetas[sample(1:nrow(thetas), size = knobj$global_parameters$final_sample_design),]
		
		if(!is.null(file_to_save)){
			saveRDS(knobj, file_to_save)
		}
		
		## Choose next experiment and generate data
		if(verbose){
			print(paste("Get data", k))
		}
		exps <- exps[exps$Cost <= credits,]
		next_exp <- sample(1:nrow(exps), size = 1) 
				
		knobj$experiments <- c(knobj$experiments, paste(exps$exp[next_exp], exps$Measurement[next_exp]) )
		
		## Simulate data, add noise and take the subset corresponding to the chosen quantity to observe
		exp_fun_next_exp <- experiment_list1[[which(names(experiment_list1) == exps[next_exp, 3])]]
		data_next_exp <- simulate_experiment(knobj$global_parameters$true_params_T, knobj, exp_fun_next_exp)
		data_next_exp <- add_noise(data_next_exp)		
		to_observe <- observables[[ as.character(exps[next_exp,1]) ]]$obs
		time_res <- observables[[ as.character(exps[next_exp,1]) ]]$reso
		knobj$datas[[length(knobj$datas) + 1]] <- list(manip =  exp_fun_next_exp, data = data_next_exp[data_next_exp[,1] %in% time_res,to_observe] )

		
		
		credits <- credits - exps[next_exp,]$Cost
		exps <- exps[-next_exp,]
		next_it <- (credits >  min(sapply(observables, FUN= function(x){x$cost})))
	
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
	risks <- c()
	if(!is.null(file_to_save)){
		saveRDS(knobj, file_to_save)
	}
	
	knobj
}
