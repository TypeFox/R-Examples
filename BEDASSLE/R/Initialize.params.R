Initialize.params <-
function(counts,sample_sizes,k,loci){
		allele.frequencies_hat <- counts/sample_sizes
			mean_allele.frequencies <- colSums(allele.frequencies_hat,na.rm=TRUE)/k
				mean_allele.frequencies[which(mean_allele.frequencies==0)] <- 0.01
				mean_allele.frequencies[which(mean_allele.frequencies==1)] <- 0.99
		mu_hat <- -log((1/mean_allele.frequencies)-1)
			sd_mu_hat <- sd(mu_hat)
			beta_hat <- 1/(sd_mu_hat)^2			
		mu_hat <- matrix(mu_hat,nrow=k,ncol=loci,byrow=TRUE)
			allele.frequencies_hat[which(allele.frequencies_hat == 0)] <- 0.01
			allele.frequencies_hat[which(allele.frequencies_hat == 1)] <- 0.99
			allele.frequencies_hat[which(is.na(allele.frequencies_hat))] <- 0.01			
		thetas_hat <- -log((1/allele.frequencies_hat)-1)-mu_hat
		initial_param_estimates <- list(thetas_hat,mu_hat,beta_hat)
			names(initial_param_estimates) <- c("thetas_hat","mu_hat","beta_hat")
		return(initial_param_estimates)
	}
