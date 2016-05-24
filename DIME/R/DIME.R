DIME <-
function(data, avg = NULL, gng.K = 2, gng.weights = NULL,
  gng.weights.cutoff = -1.345, gng.pi = NULL, 
  gng.mu = NULL, gng.sigma = NULL, gng.beta = NULL, gng.tol=1e-5, 
  gng.max.iter=2000, gng.th= NULL, gng.rep = 15, gng.fdr.cutoff = 0.1, 
  gng.sigma.diff.cutoff = NULL, gng.mu.diff.cutoff = NULL, gng.var.thres=1e2,
  gng.min.sd = NULL,
  inudge.K = 2, inudge.weights = NULL, 
  inudge.weights.cutoff = -1.345, inudge.pi = NULL, inudge.mu = NULL,
  inudge.sigma = NULL, inudge.tol=1e-5, inudge.max.iter = 2000, inudge.z=NULL,
  inudge.rep = 15, inudge.fdr.cutoff = 0.1, inudge.sigma.diff.cutoff = NULL,
  inudge.mu.diff.cutoff = NULL,  inudge.var.thres=1e2, inudge.min.sd = NULL,
  nudge.z = NULL, nudge.tol = 1e-5, nudge.max.iter = 2000, nudge.mu = NULL, 
  nudge.sigma = NULL, nudge.rep = 15, nudge.fdr.cutoff = 0.1, 
  nudge.weights = NULL,
  nudge.weights.cutoff = -1.345, nudge.pi=NULL){
  
  obs <- unlist(data);
  n <- length(obs);
 
  # checking input
  if((!is.null(gng.mu)&&is.null(dim(gng.mu)))||
      (!is.null(gng.mu)&&!is.null(dim(gng.mu))&&dim(gng.mu)[2]<gng.K)){
      return(cat('Error: Please input values for gng.mu as a matrix, 
        dim(gng.mu)[2] must be = gng.K\n'));
  }
  if((!is.null(gng.sigma)&&is.null(dim(gng.sigma)))||
      (!is.null(gng.sigma)&&!is.null(dim(gng.sigma))&&dim(gng.sigma)[2]<gng.K)){
      return(cat('Error: Please input values for gng.sigma as a matrix,
        dim(gng.sigma)[2] must be = gng.K\n'));
  }
  if((!is.null(gng.pi)&&is.null(dim(gng.pi)))||
      (!is.null(gng.pi)&&!is.null(dim(gng.pi))&&dim(gng.pi)[2]<(gng.K+2))){
      return(cat('Error: Please input values for gng.pi as a matrix, 
        dim(gng.pi)[2] must be = (gng.K+2)\n'));
  }
  if((!is.null(gng.beta)&&is.null(dim(gng.beta)))||
      (!is.null(gng.beta)&&!is.null(dim(gng.beta))&&dim(gng.beta)[2]< 2)){
      return(cat('Error: Please input values for gng.beta as a matrix,
        dim(gng.beta)[2] must be = 2\n'));
  }
  if((!is.null(gng.th)&&is.null(dim(gng.th)))||
    (!is.null(gng.th)&&!is.null(dim(gng.th))&&dim(gng.th)[2]< 2)){
      return(cat('Error: Please input values for gng.th matrix
        dim(gng.th)[2] must be = 2\n'));
  }
 if(!is.null(gng.weights)&&is.null(avg)){
      return(cat('Error: When using weights, please input the mean 
        (or log intensities) of data\n'));
 }
 if((!is.null(gng.weights)&&!is.character(gng.weights)&&
      is.null(dim(gng.weights)))||
      (!is.null(gng.weights)&&!is.null(dim(gng.weights))&&dim(gng.weights)[2]%%n!= 0)){
      return(cat('Error: Please input values for gng.weights as a matrix whose 
        no. columns = of length(unlist(data))\n'));
  }
  
  # checking input
  if((!is.null(inudge.mu)&&is.null(dim(inudge.mu)))||
    (!is.null(inudge.mu)&&!is.null(dim(inudge.mu))&&dim(inudge.mu)[2]<inudge.K)){
      return(cat('Error: Please input values for inudge.mu as a matrix, 
        dim(inudge.mu)[2] must be = inudge.K\n'));
  }
  if((!is.null(inudge.sigma)&&is.null(dim(inudge.sigma)))||
    (!is.null(inudge.sigma)&&!is.null(dim(inudge.sigma))&&dim(inudge.sigma)[2]<inudge.K)){
      return(cat('Error: Please input values for inudge.sigma as a matrix,
        dim(inudge.sigma)[2] must be = inudge.K\n'));
  }
  if((!is.null(inudge.pi)&&is.null(dim(inudge.pi)))||
    (!is.null(inudge.pi)&&!is.null(dim(inudge.pi))&&dim(inudge.pi)[2]<(inudge.K+1))){
      return(cat('Error: Please input values for inudge.pi as a matrix,
        dim(inudge.pi)[2] must be = (inudge.K+1)\n'));
  }
  if(!is.null(inudge.weights)&&is.null(avg)){
      return(cat('Error: When using weights, please input the mean 
        (or log intensities) of data\n'));
  }
  if((!is.null(inudge.weights)&&!is.character(inudge.weights)&&
      is.null(dim(inudge.weights)))||
      (!is.null(inudge.weights)&&!is.null(dim(inudge.weights))&&dim(inudge.weights)[2]%%n!= 0)){
      return(cat('Error: Please input values for inudge.weights as a matrix whose 
        no of columns = length(unlist(data))\n'));
  }
  if((!is.null(inudge.z)&&is.null(dim(inudge.z)))||
    (!is.null(inudge.z)&&!is.null(dim(inudge.z))&&(dim(inudge.z)[2]!= 2||dim(inudge.z)[1]!=n ||
    (length(which(rowSums(inudge.z)!=1))!=0)))){
      return(cat('Error: Please input values for inudge.z as a matrix with no. rows
        = length(unlist(data)) and no. columns = 2. Each row must sum to 1\n'));
  }
    
  # checking input
  if((!is.null(nudge.mu)&&is.null(dim(nudge.mu)))||
    (!is.null(nudge.mu)&&!is.null(dim(nudge.mu))&&dim(nudge.mu)[2]< 1)){
      return(cat('Error: Please input values for nudge.mu as a matrix, 
        dim(nudge.mu)[2] must be = 1\n'));
  }
  if((!is.null(nudge.sigma)&&is.null(dim(nudge.sigma)))||
    (!is.null(nudge.sigma)&&!is.null(dim(nudge.sigma))&&dim(nudge.sigma)[2]< 1)){
      return(cat('Error: Please input values for nudge.sigma as a matrix, 
        dim(nudge.sigma)[2] must be = 1\n'));
  }
  if((!is.null(nudge.pi)&&is.null(dim(nudge.pi)))||
    (!is.null(nudge.pi)&&!is.null(dim(nudge.pi))&&dim(nudge.pi)[2]<2)){
      return(cat('Error: Please input values for nudge.pi as a matrix, 
        dim(nudge.pi)[2] must be = 2\n'));
  }
  if(!is.null(nudge.weights)&&is.null(avg)){
      return(cat('Error: When using weights, please input the mean 
        (or log intensities) of data\n'));
  }
  if((!is.null(nudge.weights)&&!is.character(nudge.weights)&&
    is.null(dim(nudge.weights)))||
    (!is.null(nudge.weights)&&!is.null(dim(nudge.weights))&&dim(nudge.weights)[2]%%n!= 0)){
      return(cat('Error: Please input values for nudge.weights a matrix whose
        no. columns = length(unlist(data))\n.'));
  }
 if((!is.null(nudge.z)&&is.null(dim(nudge.z)))||
    (!is.null(nudge.z)&&!is.null(dim(nudge.z))&&(dim(nudge.z)[2]!= 2 ||dim(nudge.z)[1]!=n ||
    (length(which(rowSums(nudge.z)!=1))!=0)))){
      return(cat('Error: Please input values for nudge.z as a matrix with no. rows
        = length(unlist(data)) and no. columns = 2\n. Each row must sum to 1\n'));
  }
  
   # run NUDGE for rep times (different random number) to get best NUDGE
  bestNudgeBIC <- (-Inf);
  for (i in 1:nudge.rep){
      if((!is.null(nudge.mu) && i>dim(nudge.mu)[1])||is.null(nudge.mu))
            { mu_i <- NULL;}else{ mu_i <- nudge.mu[i,1]}
      if((!is.null(nudge.sigma) && i>dim(nudge.sigma)[1])||is.null(nudge.sigma))
            { sigma_i <- NULL;}else{ sigma_i <- abs(nudge.sigma[i,1])}
      if((!is.null(nudge.pi) && i>dim(nudge.pi)[1])||is.null(nudge.pi))
            { pi_i <- NULL;}else{ pi_i <- nudge.pi[i,1]/sum(nudge.pi[i,1])}
      if(!is.null(nudge.weights) && is.character(nudge.weights)){
      	weights <- match.arg(tolower(nudge.weights),c("lower","upper","full"));
	       weights_i <- switch(weights,
			   lower = huber(avg, nudge.weights.cutoff,'lower'),
			   upper = huber(avg, nudge.weights.cutoff,'upper'),
			   full = huber(avg, nudge.weights.cutoff,'full'))
      }
      if((!is.null(nudge.weights) && !is.character(nudge.weights)
        && i>dim(nudge.weights)[1]) ||is.null(nudge.weights)){
          weights_i <- NULL;
        }else if(!is.character(nudge.weights)){
          weights_i <- nudge.weights[i,];
        }
     if(is.null(nudge.z))
            { z_i <- NULL;}else{z_i <- nudge.z}
      nudge <- inudge.fit(obs, K = 1, weights = weights_i, 
      pi = pi_i, mu = mu_i, sigma = sigma_i, tol = nudge.tol, 
      max.iter= nudge.max.iter, z = z_i);
      # check for non-convergence
  if(is.nan(nudge$BIC) || is.infinite(nudge$BIC))	nudge$BIC <- (-Inf);
      if (bestNudgeBIC <= nudge$BIC) {
    		bestNudge <- nudge ;
  		  bestNudgeBIC <- nudge$BIC ;
  		  bestNudge$name <- "NUDGE"
     }
  }
  
  # run GNG for rep x K times (different random number) to get best GNG
  bestGngBIC <- (-Inf);
  for (i in 1:gng.rep){
  	for (k in 1:gng.K){
      if((!is.null(gng.mu) && i>dim(gng.mu)[1])||is.null(gng.mu))
            { mu_i <- NULL;}else{ mu_i <- gng.mu[i,1:k]}
      if((!is.null(gng.sigma) && i>dim(gng.sigma)[1])||is.null(gng.sigma))
            { sigma_i <- NULL;}else{ sigma_i <- abs(gng.sigma[i,1:k])}
      if((!is.null(gng.th) && i>dim(gng.th)[1])||is.null(gng.th))
            { th_i <- NULL;}else{ th_i <- gng.th[i,]}
      if((!is.null(gng.pi) && i>dim(gng.pi)[1])||is.null(gng.pi))
            { pi_i <- NULL;}else{ pi_i <- (gng.pi[i,1:(k+2)])/sum(gng.pi[i,1:(k+2)])}
      if(!is.null(gng.weights) && is.character(gng.weights)){
      	weights <- match.arg(tolower(gng.weights),c("lower","upper","full"));
		    weights_i <- switch(weights,
			   lower = huber(avg, gng.weights.cutoff,'lower'),
			   upper = huber(avg, gng.weights.cutoff,'upper'),
			   full = huber(avg, gng.weights.cutoff,'full'))
      }      
      if((!is.null(gng.weights) && !is.character(gng.weights)
            && i>dim(gng.weights)[1])||is.null(gng.weights)){
              weights_i <- NULL;
            }else if( !is.character(gng.weights)){
              weights_i <- gng.weights[i,]}
      if((!is.null(gng.beta) && i>dim(gng.beta)[1])||is.null(gng.beta))
            { beta_i <- NULL;}else{ beta_i <- gng.beta[i,1:2]}
     gng <- gng.fit(obs, K = k, tol=gng.tol, max.iter=gng.max.iter,
          mu= mu_i, sigma=sigma_i, pi=pi_i, th= th_i, beta=beta_i,
          weights = weights_i);
      if(is.null(gng.min.sd)) gng.min.sd = 0.1*(sd(obs));
  		# check for AIC = NaN
  		if(is.nan(gng$BIC) || is.infinite(gng$BIC)||
        (max(gng$sigma)^2/min(gng$sigma)^2) > gng.var.thres ||
        min(gng$sigma) < gng.min.sd)	gng$BIC <- (-Inf);
  		if (bestGngBIC <= gng$BIC) {
  			bestGng <- gng ;
  			bestGngBIC <- gng$BIC;
				bestGng$name <- "GNG";
  		}
  	}
  }
  
  # run iNudge for rep x K times (different random number) to get best iNUDGE
  bestiNudgeBIC <- (-Inf);
  for (i in 1:inudge.rep){
  	for (k in 1:inudge.K){
      if((!is.null(inudge.mu) && i>dim(inudge.mu)[1])||is.null(inudge.mu))
            { mu_i <- NULL;}else{ mu_i <- inudge.mu[i,1:k]}
      if((!is.null(inudge.sigma) && i>dim(inudge.sigma)[1])||is.null(inudge.sigma))
            { sigma_i <- NULL;}else{ sigma_i <- abs(inudge.sigma[i,1:k])}
      if((!is.null(inudge.pi) && i>dim(inudge.pi)[1])||is.null(inudge.pi))
            { pi_i <- NULL;}else{ pi_i <- (inudge.pi[i,1:(k+1)])/sum(inudge.pi[i,1:(k+1)])}
      if(!is.null(inudge.weights) && is.character(inudge.weights)){
      	 weights <- match.arg(tolower(inudge.weights),c("lower","upper","full"));
		     weights_i <- switch(weights,
			   lower = huber(avg, inudge.weights.cutoff,'lower'),
			   upper = huber(avg, inudge.weights.cutoff,'upper'),
			   full = huber(avg, inudge.weights.cutoff,'full'))
      }
      if((!is.null(inudge.weights) && !is.character(inudge.weights) 
            && i>dim(inudge.weights)[1])||is.null(inudge.weights)){
              weights_i <- NULL;
            }else if(!is.character(inudge.weights)){
              weights_i <- inudge.weights[i,]
            }
      if(is.null(inudge.z))
            {z_i <- NULL;}else{ z_i <- inudge.z}
     inudge <- inudge.fit(obs, K = k,  weights = weights_i, pi = pi_i, mu= mu_i,
            sigma=sigma_i, tol=inudge.tol, max.iter=inudge.max.iter, z=z_i);
      if(is.null(inudge.min.sd)) inudge.min.sd = 0.1*(sd(obs));
  		# check for AIC = NaN
  		if(is.nan(inudge$BIC) || is.infinite(inudge$BIC) ||
      (max(inudge$sigma)^2/min(inudge$sigma)^2) > inudge.var.thres ||
      min(inudge$sigma) < inudge.min.sd)	inudge$BIC <- (-Inf);
  		if (bestiNudgeBIC <= inudge$BIC) {
  			bestiNudge <- inudge ;
  			bestiNudgeBIC <- inudge$BIC;
  			bestiNudge$name <- "iNUDGE";
  		}
  	}   
  }
  # skip classification if BIC = -Inf
  if(bestNudge$BIC!=-Inf){
  bestNudge <- DIME.classify(data,bestNudge,nudge.fdr.cutoff);
  }else{ cat('Warning! NUDGE model did not converge\n'); bestNudge$class = rep(0,n);}
  if(bestGng$BIC!=-Inf){
  bestGng <- DIME.classify(data,bestGng,gng.fdr.cutoff,gng.sigma.diff.cutoff,
    gng.mu.diff.cutoff);
    }else{cat('Warning! GNG model did not converge\n');bestGng$class = rep(0,n);}
  if(bestiNudge$BIC!=-Inf){  
  bestiNudge <- DIME.classify(data,bestiNudge,inudge.fdr.cutoff,
    inudge.sigma.diff.cutoff, inudge.mu.diff.cutoff);
    }else{cat('Warning! iNUDGE model did not converge\n');bestiNudge$class = rep(0,n);}
  
  # model selection
  if(max(bestiNudge$AIC,bestGng$AIC,bestNudge$AIC)==bestGng$AIC){
      bestMixture <- bestGng;
  }else if(max(bestiNudge$AIC,bestGng$AIC,bestNudge$AIC)==bestiNudge$AIC){
      bestMixture <- bestiNudge;
  }else{
      bestMixture <- bestNudge;
  }
  mixture <- list(best=bestMixture,nudge=bestNudge,gng=bestGng,inudge=bestiNudge);
  return (mixture);
}

