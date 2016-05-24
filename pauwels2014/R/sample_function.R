sample_function <-
function(knobj){

	# Initialization
	thetasT <- c()

	# Repeat for each mode
	for(i in 1:knobj$global_parameters$n_multi_mod_weight){
		FAIL <- TRUE
		print("optim")

		# Consider "good" samples, which give correct fit and no NA
		while(FAIL){

			# Initialize
			theta <- runif( min = 0, max = 50, length(knobj$global_parameters$param_names) ) + 25
			names(theta) <- knobj$global_parameters$param_names

			# Find a mode
			temp <- BFGS_special(theta, knobj, eval_log_like_knobj, verbose = F)
			theta <- temp$theta

			# Filter for parameters with bad fit (heuristic tests)
			FAIL <- temp$fail | (eval_log_like_knobj(theta, knobj) < knobj$global_parameters$max_log_like) # good likelihood

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

			)) > knobj$global_parameters$centrality_ratio # corresponding curves should stand "in the middle" of the data
			FAIL <- FAIL | tempF
		}

		# Sample starting drom this mode a,d take a subsample
		print("sample")
		thetas <- generate_sample(theta, knobj, N = (knobj$global_parameters$sample_burn_in + knobj$global_parameters$sample_to_keep1), step = 1, verbose = F)
		thetas <- thetas[-(1:knobj$global_parameters$sample_burn_in),][sample(1:knobj$global_parameters$sample_to_keep1, size = knobj$global_parameters$sample_to_keep1),]
		thetasT[[length(thetasT)+1]] <- thetas
	}

	# Reoptimize from previous samples
	thetas_prev <- c()
	if(length(knobj$datas) > 1){
		thetas_prev <- knobj$datas[[length(knobj$datas)-1]]$thetas
	}

	if(!is.null(thetas_prev)){
		to_ressample <- sample(1:nrow(thetas_prev))
		for(i in 1:knobj$global_parameters$n_multi_mod_weight){

			print("re-optim")

			# Initialize
			theta <- thetas_prev[to_ressample[i],]
			names(theta) <- knobj$global_parameters$param_names

			# Find a mode
			temp <- BFGS_special(theta, knobj, eval_log_like_knobj, verbose = F)
			theta <- temp$theta

			# Filter for parameters with bad fit (heuristic tests)
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

			# Sample from this mode
			if(!FAIL){
				print("sample")
				thetas <- generate_sample(theta, knobj, N = (knobj$global_parameters$sample_burn_in + knobj$global_parameters$sample_to_keep1), step = 1, verbose = F)
				thetas <- thetas[-(1:knobj$global_parameters$sample_burn_in),][sample(1:knobj$global_parameters$sample_to_keep1, size = knobj$global_parameters$sample_to_keep1),]
				thetasT[[length(thetasT)+1]] <- thetas
			}
		}
	}
	#lls <- apply(thetasT, 1, eval_log_like_knobj, knobj = knobj, fit = T)
	n_mode <- length(thetasT)
	thetasT <- sapply(1:n_mode, function(k){thetasT[[k]]}, simplify ="array")
	
	mu_est <- t(sapply(1:n_mode, function(k){
		apply(thetasT[,,k], 2, mean)
	}))

	sigma_est <- sapply(1:n_mode, function(k){
		var(thetasT[,,k])
	}, simplify = "array")
	
	lls <- apply(mu_est, 1, eval_log_like_knobj, knobj = knobj, fit = T)
	lls <- exp(lls + min(lls))
	lls <- lls / mean(lls)
	
	matA <- t(sapply(1:n_mode, function(i){  
		sapply(1:n_mode, function(j){
				dmvnorm(mu_est[i,], mu_est[j,], sigma_est[,,j])
		})
	}))
	print(matA)
	matA <- matA / mean(matA)
	
	pi_est <- solve(matA + diag(1e-10,nrow(matA))) %*% lls
	pi_est <- pi_est / sum(pi_est)

	pi_est2 <- proj_grad(matA, lls, pi_est)
	pi_est2 <- pi_est2 / sum(pi_est2)
	if(sum(is.na(pi_est2)) > 0){
		pi_est2 <- rep(1/length(pi_est2), length(pi_est2))
	}
	
	##lls <- log_normalize(lls)[[1]]
	thetasT <- thetasT[sample(1:nrow(thetasT), size = knobj$global_parameters$final_sample, replace = T), ,]
	
	temp_select <- sample(1:n_mode, size = knobj$global_parameters$final_sample, replace = T, prob = pi_est2)
	
	thetasT <- t(sapply(1:length(temp_select), function(k){thetasT[k,,temp_select[k]]}))
	thetasT
}
