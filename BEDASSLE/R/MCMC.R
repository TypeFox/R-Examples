MCMC <-
function(counts,
				sample_sizes,
				D,
				E,
				k,
				loci,
				delta,
				aD_stp,
				aE_stp,
				a2_stp,
				thetas_stp,
				mu_stp,
				ngen,
				printfreq,
				savefreq,
				samplefreq,
				directory=NULL,
                prefix="",
				continue=FALSE,
				continuing.params=NULL){
	
	if(!is.null(directory)){
		setwd(directory)
	}	
	if(min(D[upper.tri(D)]) < delta){
		warning("the value chosen for delta may too large compared to the geographic distances between the sampled populations")
	}
    if (is.matrix(E)) { E <- list(E) }  # allow E to come not in a list if there's only one
	if(any(lapply(E,function(EE) min(EE[upper.tri(EE)])) < delta)){
		warning("the value chosen for delta may too large compared to the ecological distances between the sampled populations")
	}
	
	#Declare variables
		LnL_thetas <- numeric(ngen/samplefreq) 				
		LnL_thetas_vec <- numeric(loci)						
		LnL_counts <- numeric(ngen/samplefreq)				
		LnL <- numeric(ngen/samplefreq)						
		Prob <- numeric(ngen/samplefreq)					
		a0 <- numeric(ngen/samplefreq)						
		aD <- numeric(ngen/samplefreq)						
		aE <- matrix(0,nrow=length(E),ncol=ngen/samplefreq)	
		a2 <- numeric(ngen/samplefreq)						
		beta <- numeric(ngen/samplefreq)					
		a0_moves <- numeric(ngen/samplefreq)				
		aD_moves <- numeric(ngen/samplefreq)				
		aE_moves <- numeric(ngen/samplefreq)				
		a2_moves <- numeric(ngen/samplefreq)				
		thetas_moves <- numeric(ngen/samplefreq)			
		mu_moves <- numeric(ngen/samplefreq)				
		beta_moves <- numeric(ngen/samplefreq)				
		aD_accept <- numeric(ngen/samplefreq)				
		aE_accept <- numeric(ngen/samplefreq)				
		a2_accept <- numeric(ngen/samplefreq)				
		thetas_accept <- numeric(ngen/samplefreq)			
		mu_accept <- numeric(ngen/samplefreq)				
		
	if(!continue) {
		#INITIALIZE MCMC
				Prob[1] <- -Inf
				covariance <- matrix(0,nrow=k,ncol=k)
				pos.def.counter <- 0

			while(Prob[1] == -Inf){
				while(!is.positive.definite(covariance) && pos.def.counter < 100){
						a0[1] <- runif(1,0,4)												
						aD[1] <- runif(1,0,4)												
						aE[,1] <- runif(length(E),0,4)										
						a2[1] <- runif(1,0,2)												
					covariance <- Covariance(a0[1],aD[1],aE[,1],a2[1],D,E,delta)	
					initial.params <- Initialize.params(counts,sample_sizes,k,loci) 		
					beta[1] <- initial.params$beta_hat
					mu <- initial.params$mu_hat
					thetas <- initial.params$thetas_hat
						allele.frequencies <- transform_frequencies(thetas,mu)
					Initial.parameters <- list(a0=a0[1],aD=aD[1],aE=aE[,1],a2=a2[1],beta=beta[1],mu=mu,thetas=thetas)
					pos.def.counter <- pos.def.counter + 1			
				}
					if(pos.def.counter > 99){
						stop("the initial covariance matrix is not positive definite! Either attempt to re-initialize MCMC, or increase the size of the delta shift")
					}			
                save(Initial.parameters,file=paste(prefix,"Initial.parameters.Robj",sep=''))						
				LnL_thetas_vec <- Likelihood_thetas(thetas,covariance)								
				LnL_thetas[1] <- sum(Likelihood_thetas(thetas,covariance))
				LnL_counts_mat <- Likelihood_counts(counts,sample_sizes,allele.frequencies)
				LnL_counts[1] <- sum(LnL_counts_mat)
				LnL[1] <- LnL_thetas[1] + LnL_counts[1]
				prior_prob_beta <- Prior_prob_beta(beta[1])
				prior_prob_mu <- Prior_prob_mu(mu,beta[1])
				prior_prob_alpha0 <- Prior_prob_alpha0(a0[1])
				prior_prob_alphaD <- Prior_prob_alphaD(aD[1])
				prior_prob_alphaE <- Prior_prob_alphaE(aE[,1])
				prior_prob_alpha2 <- Prior_prob_alpha2(a2[1])				
				Prob[1] <- LnL[1] + prior_prob_alpha0 + prior_prob_alphaD + prior_prob_alphaE + prior_prob_alpha2 + sum(prior_prob_mu) + prior_prob_beta
			}
	} else {
		#Choose Initial Parameter Values
			load(continuing.params)
			a0[1] <- continuing.params$a0											
			aD[1] <- continuing.params$aD											
			aE[,1] <- continuing.params$aE											
			a2[1] <- continuing.params$a2											
				covariance <- Covariance(a0[1],aD[1],aE[,1],a2[1],D,E,delta)	
				if(!is.positive.definite(covariance)){ 								
					stop("the initial covariance matrix is not positive definite! Either attempt to re-initialize MCMC, or increase the size of the delta shift")
				}			
			beta[1] <-  continuing.params$beta											
			mu <-  continuing.params$mu													
			thetas <-  continuing.params$thetas											
			allele.frequencies <- transform_frequencies(thetas,mu)
                Initial.parameters <- list(a0=a0[1],aD=aD[1],aE=aE[,1],a2=a2[1],beta=beta[1],mu=mu,thetas=thetas)
			save(Initial.parameters,file=paste(prefix,"Initial.parameters.Robj",sep=''))				
	}
	
	#Calculate Initial Likelihood
		LnL_thetas_vec <- Likelihood_thetas(thetas,covariance)								
		LnL_thetas[1] <- sum(Likelihood_thetas(thetas,covariance))
		LnL_counts_mat <- Likelihood_counts(counts,sample_sizes,allele.frequencies)
		LnL_counts[1] <- sum(LnL_counts_mat)
		LnL[1] <- LnL_thetas[1] + LnL_counts[1]
		prior_prob_beta <- Prior_prob_beta(beta[1])
		prior_prob_mu <- Prior_prob_mu(mu,beta[1])
		prior_prob_alpha0 <- Prior_prob_alpha0(a0[1])
		prior_prob_alphaD <- Prior_prob_alphaD(aD[1])
		prior_prob_alphaE <- Prior_prob_alphaE(aE[1])
		prior_prob_alpha2 <- Prior_prob_alpha2(a2[1])				
		Prob[1] <- LnL[1] + prior_prob_alpha0 + prior_prob_alphaD + prior_prob_alphaE + prior_prob_alpha2 + sum(prior_prob_mu) + prior_prob_beta
			if(!is.finite(Prob[1])){													
				stop("Initial probability of model is NEGATIVE INFINITY! Please attempt to initiate chain again.")
			}
	last.params <- list(a0[1],aD[1],aE[,1],a2[1],beta[1],delta,covariance,thetas,mu,allele.frequencies,aD_stp,aE_stp,a2_stp,thetas_stp,mu_stp,k,LnL_thetas_vec,LnL_counts_mat,prior_prob_mu,prior_prob_alpha0,prior_prob_alphaD,prior_prob_alphaE,prior_prob_alpha2,prior_prob_beta,a0_moves[1],aD_moves[1],aE_moves[1],a2_moves[1],thetas_moves[1],mu_moves[1],beta_moves[1],aD_accept[1],aE_accept[1],a2_accept[1],thetas_accept[1],mu_accept[1],counts,sample_sizes,loci,D,E)
	names(last.params) <- c("a0","aD","aE","a2","beta","delta","covariance","thetas","mu","allele.frequencies","aD_stp","aE_stp","a2_stp","thetas_stp","mu_stp","k","LnL_thetas_vec","LnL_counts_mat","prior_prob_mu","prior_prob_alpha0","prior_prob_alphaD","prior_prob_alphaE","prior_prob_alpha2","prior_prob_beta","a0_moves","aD_moves","aE_moves","a2_moves","thetas_moves","mu_moves","beta_moves","aD_accept","aE_accept","a2_accept","thetas_accept","mu_accept","counts","sample_sizes","loci","D","E")

	#Run the MCMC
	Updates <- list(Update_a0,Update_aD,Update_aE,Update_a2,Update_thetas,Update_mu,Update_beta)
	
	for(i in 2:ngen) {
		x <- sample(c(1:length(Updates)),1)
		new.params <- Updates[[x]](last.params)
        if (length(new.params$aE)<length(E)) { browser() }

		if(i%%samplefreq == 0){
			j <- i/samplefreq
			a0[j] <- new.params$a0
			aD[j] <- new.params$aD
			aE[,j] <- new.params$aE	
			a2[j] <- new.params$a2
			beta[j] <- new.params$beta
			LnL_thetas[j] <- sum(new.params$LnL_thetas_vec)
			LnL_counts[j] <- sum(new.params$LnL_counts_mat)
			LnL[j] <- LnL_thetas[j]+LnL_counts[j]
			Prob[j] <- LnL[j]+
						sum(new.params$prior_prob_mu)+
						new.params$prior_prob_beta+
						new.params$prior_prob_alpha0+
						new.params$prior_prob_alphaD+
						new.params$prior_prob_alphaE+
						new.params$prior_prob_alpha2
			a0_moves[j] <- new.params$a0_moves
			aD_moves[j] <- new.params$aD_moves
			aE_moves[j] <- new.params$aE_moves
			a2_moves[j] <- new.params$a2_moves
			thetas_moves[j] <- new.params$thetas_moves
			mu_moves[j] <- new.params$mu_moves
			beta_moves[j] <- new.params$beta_moves		
			aD_accept[j] <- new.params$aD_accept
			aE_accept[j] <- new.params$aE_accept
			a2_accept[j] <- new.params$a2_accept
			thetas_accept[j] <- new.params$thetas_accept
			mu_accept[j] <- new.params$mu_accept
		}
		
		last.params <- new.params

		if(i%%printfreq == 0){
			P <- (sum(new.params$LnL_thetas_vec)
					+ sum(new.params$LnL_counts_mat)
					+ sum(new.params$prior_prob_mu)
					+ new.params$prior_prob_beta
					+ new.params$prior_prob_alpha0
					+ new.params$prior_prob_alphaD
					+ new.params$prior_prob_alphaE
					+ new.params$prior_prob_alpha2)
			print(i)
			print(P)
		}
				
		if(i%%savefreq == 0){	
			save(last.params,
				LnL_thetas,
				LnL_counts,
				LnL,Prob,
				a0,aD,aE,a2,beta,
				samplefreq,ngen,
				a0_moves,aD_moves,aE_moves,a2_moves,thetas_moves,mu_moves,beta_moves,
				aD_accept,aE_accept,a2_accept,thetas_accept,mu_accept,
				aD_stp,aE_stp,a2_stp,thetas_stp,mu_stp,
				file=paste(prefix,sprintf("MCMC_output%d.Robj",1),sep=''))
		}
	}
    return(paste("Output",i,"runs to",paste(prefix,"MCMC_output*.Robj",sep=''),"."))
}
