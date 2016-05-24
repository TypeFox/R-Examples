sample_function_single_mod <-
function(knobj){
	
	thetasT <- c()
	FAIL <- TRUE
	print("optim")
	while(FAIL){
    	# Initialize
		theta <- runif( min = 0, max = 50, length(knobj$global_parameters$param_names) ) + 25
		names(theta) <- knobj$global_parameters$param_names
		
		# Find a mode
		temp <- BFGS_special(theta, knobj, eval_log_like_knobj)
		theta <- temp$theta
		
		# Filter for good samples
		FAIL <- temp$fail | (eval_log_like_knobj(theta, knobj) < knobj$global_parameters$max_log_like)
		tempF <- max(abs(
			sapply( 1:length(knobj$datas), 
				FUN = function(i){
					if(!is.null(knobj$datas[[i]]$data)){
						temp <- simulate_experiment(theta, knobj, knobj$datas[[i]]$manip);
						mean( sign( ( knobj$datas[[i]]$data - temp[temp[,1] %in% knobj$datas[[i]]$data[,1], dimnames(temp)[[2]] %in% dimnames(knobj$datas[[i]]$data)[[2]]]
						)[,-1]))
					}else{
						0
					}
				}
			)

		)) > knobj$global_parameters$centrality_ratio
		FAIL <- FAIL | tempF
	}
	# Generate a sample initialized at this mode
	print("sample")
	thetas <- generate_sample(theta, knobj, N = (knobj$global_parameters$sample_burn_in + knobj$global_parameters$sample_to_keep1), step = 1)
	thetas <- thetas[-(1:knobj$global_parameters$sample_burn_in),][sample(1:knobj$global_parameters$final_sample, size = knobj$global_parameters$final_sample),]
}
