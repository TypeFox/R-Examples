ABC_SMC_nLTT <- function(tree,statistics,simFunc,initEpsilon,PRIOR_GEN,PRIOR_DENS,numParticles,sigma=0.05,stopRate=1e-5) {

 parameters <- PRIOR_GEN; #just to get the number of parameters to be estimated.
 # compute the observed statistics
 obs_statistics <- c()
  for (i in 1:length(statistics)) {
    obs_statistics[i] <- statistics[[i]](tree)
  }
  
  stats <- c()
  
  
  #generate a matrix with epsilon values, we assume that the SMC algorithm converges within 100 iterations
  epsilon <- matrix(nrow=50,ncol=length(initEpsilon));
  for(j in 1:length(initEpsilon)) {
   for(i in 1:50) {
    epsilon[i,j] <- initEpsilon[j] * exp(-0.5 * (i-1));
   }
  }
  
  #store weights
  newWeights <- c();
  newParams <- list(c(1:length(parameters)));
  previousWeights <- c();
  previousParams  <- list(c(1:length(parameters)));
  indices <- 1:numParticles;

  for( i in 1:50 ) { #convergence is expected within 100 iterations, usually convergence occurs within 10 iterations
    cat("\nGenerating Particles for iteration\t",i,"\n")
    cat("0--------25--------50--------75--------100\n")
    cat("*"); flush.console();
	
    PRINT_FREQ <- 20;
    tried <- 0;
    numberAccepted <- 0;


    if(i > 1) { #replace all vectors
     previousWeights <- newWeights / sum(newWeights); #normalize the weights and store them as previous weights.\
     newWeights <- c(); #remove all currently stored weights
     previousParams <- newParams; #store found params
     newParams <- list(c(1:length(parameters))); #clear new params
    }
       
    while(numberAccepted < numParticles) {
	   if(i == 1) { #in this initial step, generate parameters from the prior
		  parameters <- PRIOR_GEN();
      } else { #if not in the initial step, generate parameters from the weighted previous distribution:
        
        index <- sample(x=indices,size=1,replace=TRUE,prob=previousWeights);
		for(p_index in 1:length(parameters)) {	
			parameters[p_index] <- previousParams[[index]][p_index];
		}

        toChange <- sample(1:length(parameters),1);
        eta <- log(parameters[toChange]) + rnorm(1,0,sigma); # perturb the parameter a little bit, on log scale, so parameter doesn't go < 0.
        parameters[toChange] <- exp(eta) 
     }  
  
     
     if(PRIOR_DENS(parameters) > 0) 
     {
        new_tree <- simFunc(parameters); 
        accept <- TRUE;
		   
        for(k in 1:length(statistics)) {
          stats[k] <- statistics[[k]](new_tree)
          if(is.na(stats[k])) stats[k] <- Inf;
        }
		  
        for (k in 1:length(statistics)) {
          if ( abs(stats[k]-obs_statistics[k]) > epsilon[i,k] ) {
            accept <- FALSE
            if(i == 1) accept <- TRUE  #the first step always accepts
            break
          }
        }

        if ( accept ) {
          numberAccepted <- numberAccepted + 1
          newParams[[numberAccepted]] <- parameters;
          weight_accepted = 1;
          if( i > 1) weight_accepted = calculateWeight(previousWeights,previousParams,parameters,sigma,PRIOR_DENS);
          newWeights[numberAccepted] = weight_accepted;
          
          if ((numberAccepted)%%(numParticles/PRINT_FREQ) == 0) {
                        cat("**"); flush.console(); 
          }
		}
	 }
	    
      tried <- tried + 1; #convergence if the acceptance rate gets too low
      if(tried > (1/stopRate)) { 
        if((numberAccepted / tried) < stopRate) {
			output <- c();
			for(k in 1:length(parameters)) {
			  add <- c();
			  for(m in 1:length(previousParams[[k]])) {
					add <- c(add,previousParams[[k]][m]);
			  }
			  output <- rbind(output,add);
			}
			return( output);	
			
        }
      }
    }
  }
  
  
  output <- c();
  for(k in 1:length(parameters)) {
    add <- c();
	for(m in 1:length(previousParams[[k]])) {
		add <- c(add,previousParams[[k]][m]);
	}
	output <- rbind(output,add);
  }
  return( output);	
}

calculateWeight <- function(W,P,current,S,PRIOR_DENS) {
	sum <- 0;
	vals <- c();
	for(i in 1:length(P))
    { 
      diff1 <- log(current[1]) - log(P[[i]][1]);
	    diff2 <- log(current[2]) - log(P[[i]][2]);
      vals[i] <- W[i] * dnorm(diff1,mean=0,sd=S) * dnorm(diff2,mean=0,sd=S);   
    }
   
   numerator <- PRIOR_DENS(current);
    
   return(numerator / sum(vals));
}

MCMC_nLTT <- function(phy,likelihoodFunction,parameters,logTransforms,iterations,burnin=round(iterations/3),thinning=1,sigma=1) {

  # create a list for the samples
  chain = array(dim = c(floor(iterations/thinning)+1,length(parameters))) #reserve memory for the chain, for large chains we might consider writing to a file instead of storing in memory
   

  # pre-compute current posterior probability
  pp <- likelihoodFunction(parameters,phy)

  cat("\nGenerating Particles for iteration\t",i,"\n")
  cat("0--------25--------50--------75--------100\n")
  cat("*"); flush.console();
  PRINT_FREQ <- 20;


  for (i in 1:(burnin+iterations)) {

    # propose new values
    for ( j in 1:length(parameters) ) {
      if ( logTransforms[j] == TRUE ) {
        if (parameters[j] == 0) {
          stop("Cannot propose new value for a parameter with value 0.0.")
        }
        eta           <- log(parameters[j]) ### propose a new value for parameter[j]
        new_eta       <- eta + rnorm(1,0,sigma)
        new_val       <- exp(new_eta)
        hr            <- log(new_val / parameters[j]) # calculate the Hastings ratio
        parameters[j] <- new_val
        new_pp        <- likelihoodFunction(parameters,phy)
        # accept / reject
        if ( is.finite(new_pp) && is.finite(hr) &&  new_pp-pp+hr > log(runif(1,0,1)) ) {
          pp <- new_pp
        } else {
          parameters[j] <- exp(eta)
        }
      } else {
        eta           <- parameters[j] ### propose a new value for parameter[j]
        new_val       <- eta + rnorm(1,0,sigma)
        hr            <- 0.0 # calculate the Hastings ratio
        parameters[j] <- new_val
        new_pp        <- likelihoodFunction(parameters,phy)
        # accept / reject
        if ( is.finite(new_pp) && is.finite(hr) &&  new_pp-pp+hr > log(runif(1,0,1)) ) {
          pp <- new_pp
        } else {
          parameters[j] <- eta
        }
      }

    }

    # sample the parameter
    if (i >= burnin) {
	if ((i)%%((iterations-burnin)/PRINT_FREQ) == 0) {
                        cat("**"); flush.console(); 
      }
      if ( (i-burnin) %% thinning == 0 ) {
        chain[(i-burnin)/thinning+1,] <- parameters
      }
    }
  }
  cat("Finished MCMC.\n")

  return(as.mcmc(chain)) #return a mcmc object, used by coda to plot
}















    
    
