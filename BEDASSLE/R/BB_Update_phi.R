BB_Update_phi <-
function(last.params){
		phi_prime <- last.params$phi + rnorm(last.params$k,0,last.params$phi_stp)
		prior_prob_phi_prime <- BB_Prior_prob_phi(phi_prime)
		LnL_counts_mat_prime <- BB_Likelihood_counts(phi_prime,last.params$counts,last.params$sample_sizes,last.params$allele.frequencies)
		new.params <- last.params
			for(i in 1:last.params$k){
				if(prior_prob_phi_prime[i] != -Inf){
					if(!is.na(sum(LnL_counts_mat_prime[,i]))){						
						if(exp((sum(LnL_counts_mat_prime[i,])+prior_prob_phi_prime[i])-(sum(last.params$LnL_counts_mat[i,])+last.params$prior_prob_phi[i])) >= runif(1)){
							new.params$phi[i] <- phi_prime[i]
							new.params$LnL_counts_mat[i,] <- LnL_counts_mat_prime[i,]
							new.params$prior_prob_phi[i] <- prior_prob_phi_prime[i]
							new.params$phi_accept <- new.params$phi_accept+1
						}
					}
				}	
			}	
		new.params$phi_moves <- new.params$phi_moves+new.params$k
		return(new.params)
	}
