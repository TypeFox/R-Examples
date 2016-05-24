Update_mu <-
function(last.params){
		mu_update <- rnorm(last.params$loci,0,last.params$mu_stp)
		mu_prime <- Shift(last.params$mu,mu_update)
		prior_prob_mu_prime <- Prior_prob_mu(mu_prime,last.params$beta)
		allele.frequencies_prime <- transform_frequencies(last.params$thetas,mu_prime)
		LnL_counts_mat_prime <- Likelihood_counts(last.params$counts,last.params$sample_sizes,allele.frequencies_prime)
		new.params <- last.params
			for(i in 1:ncol(last.params$mu)){
				if(exp((prior_prob_mu_prime[i]+sum(LnL_counts_mat_prime[,i])) - (last.params$prior_prob_mu[i]+sum(last.params$LnL_counts_mat[,i]))) >= runif(1)){
					new.params$mu[,i] <- mu_prime[,i]
					new.params$allele.frequencies[,i] <- allele.frequencies_prime[,i]
					new.params$LnL_counts_mat[,i] <- LnL_counts_mat_prime[,i]
					new.params$prior_prob_mu[i] <- prior_prob_mu_prime[i]
					new.params$mu_accept <- new.params$mu_accept + 1
				}
			}
		new.params$mu_moves <- new.params$mu_moves + last.params$loci
		return(new.params)
	}
