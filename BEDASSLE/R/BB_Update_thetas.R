BB_Update_thetas <-
function(last.params){
		new.params <- last.params
		cholcovmat <- chol(last.params$covariance)
		thetas_prime <- sapply(1:ncol(last.params$thetas),function(i){last.params$thetas[,i]+rnorm(last.params$k,0,last.params$thetas_stp)%*%cholcovmat})		
		allele.frequencies_prime <- transform_frequencies(thetas_prime,last.params$mu)
		LnL_thetas_vec_prime <- Likelihood_thetas(thetas_prime,last.params$covariance)
		LnL_counts_mat_prime <- BB_Likelihood_counts(last.params$phi,last.params$counts,last.params$sample_sizes,allele.frequencies_prime)
			for(i in 1:last.params$loci){
				if(!is.na(sum(LnL_counts_mat_prime[,i])) && !is.na(LnL_thetas_vec_prime[i])){
					if(exp((LnL_thetas_vec_prime[i]+sum(LnL_counts_mat_prime[,i])) - (last.params$LnL_thetas_vec[i]+sum(last.params$LnL_counts_mat[,i]))) >= runif(1)){
						new.params$thetas[,i] <- thetas_prime[,i]
						new.params$allele.frequencies[,i] <- allele.frequencies_prime[,i]
						new.params$LnL_thetas_vec[i] <- LnL_thetas_vec_prime[i]
						new.params$LnL_counts_mat[,i] <- LnL_counts_mat_prime[,i]
						new.params$thetas_accept <- new.params$thetas_accept + 1
					}
				}
			}
			new.params$thetas_moves <- new.params$thetas_moves + last.params$loci	
			return(new.params)
	}
