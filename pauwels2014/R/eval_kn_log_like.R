eval_kn_log_like <-
function(theta, initial_conditions, data, knobj, fail_incoming = F, simu = NULL, fit = F){
  ## Evaluate the log likelihood of a single experiment
  ## Log likelihood includes model terms + additional terms
  ## promoting smoothness, reasonable values etc
  
  if(is.null(simu)){
  	temp <- simulate_experiment_no_transform(theta, initial_conditions, knobj)
	}else{
		temp <- simu
	}
	
	fail <- (attr(temp,"istate")[1] == -1)
	if(fail){
		#print("Fail...")
	}
	
	temp_val <- temp[temp[,1] %in% data[,1], colnames(temp) %in% colnames(data)]
	
	res <- log_likelihood(temp, temp_val, data, fit = fit)
	
	if(fail_incoming){
		tmp <- res
		res <- c()
		res$res <- tmp
		res$fail <- fail
		res
	}
	else{
		res
	}
}
