####################################################################
######     Bayesian change point detection using metadata    #######
###### Given a model, MCMC drawls of phi, alpha, sigmasq, mu #######
####################################################################

######## Inputs #######
## X: observations, assumed to be normal with mean shifts
## meta: 1 if meta point, 0 otherwse. same length as X. meta[1] = 0
## trend: TRUE or FALSE, whether to include linear trend with time
## parameters: mu0, nu0, phi.lower, phi.upper, start.phi, sd.xi

######## Outputs ########
## Phi: 	a scalar, if EB == TRUE, 
##			(iter / thin)-length vector, if EB == FALSE
##         	MCMC samples of phi
## Sigmasq: a scalar, if EB == TRUE, 
##			(iter / thin)-length vector, if EB == FALSE
##         	MCMC samples of sigmasq
## Alpha: 	(iter / thin)-length vector, if trend == TRUE
##         	MCMC samples of sigmasq
## Mu: 		(iter / thin) * N matrix
##         	MCMC samples of Mu (length N)
## change.phi:	ratio of accepting a new phi in the MCMC chain, if EB == FALSE


####### NOTE #######
## Default strategy: use mean(X) as mu0 if not speficied 

bcpmeta.parameters = function(X, meta, eta, iter = 1e4, thin = 1e1, trend = TRUE, 
					EB = TRUE, mu0 = NULL, nu0 = 5, phi.lower = -0.99,
					phi.upper = 0.99, sd.xi = 0.1, start.phi = NULL, 
					burnin = 0.2, track.time = TRUE, show.summary = TRUE,
					start.year = 1, meta.year = FALSE, eta.year = FALSE){

  t.start = proc.time();

  ## pre-process X
  n = length(X);  
  if(length(mu0) == 0)
    mu0 = mean(X);
  Xtilde = X - mu0;
  
  ## pre-process meta
  if(length(meta) < n){
  	if(meta.year == FALSE)
  	  meta = loc2dirac(meta, n); ## time 1 cannot be a metadata point
  	if(meta.year == TRUE)
  	  meta = loc2dirac(meta - start.year + 1, n); ## time 1 cannot be a metadata point
  }
  if(length(meta) > n)
    stop('Length of metadata cannot exceed length of time series.')

  ## pre-process eta
  if(length(eta) < n){
  	if(eta.year == FALSE)
  	  eta = loc2dirac(eta, n); ## time 1 cannot be a changepoint
  	if(eta.year == TRUE)
  	  eta = loc2dirac(eta - start.year + 1, n); ## time 1 cannot be a changepoint
  }
  if(length(eta) > n)
    stop('Length of changepoint configuration cannot exceed length of time series.')
  
  ######## EB ########
  if(EB == TRUE){
  	lml.eb = lMarLik.EB.normAR(Xtilde, eta, nu0, phi.lower, phi.upper, trend);
  	Phi = lml.eb$phi.eb;
  	para.eb = up.parameters.EB.normAR(Xtilde, eta, Phi, iter = (iter / thin), nu0, trend);
  	Sigmasq = para.eb$sigmasq;
  	Mu = para.eb$Mu + mu0;
  	if(trend == TRUE)
  	  Alpha = para.eb$Alpha;
    cat(paste('100% completed.\n', sep = ''));
    cat('\n');
    change.phi = NULL;
  }
  	  
  ######## FB ########
  ## output matrices
  if(EB == FALSE){
    Phi = Sigmasq = rep(NA, nrow = round(iter / thin));  
    Mu = matrix(NA, ncol = sum(eta) + 1, nrow = round(iter / thin));
    #colnames(Mu) = paste('T', 1:n, sep = '');
    if(trend == TRUE)
 	  Alpha = rep(NA, nrow = round(iter / thin));  

    ## initial values   
    if(length(start.phi) == 0)
      phi = runif(1, phi.lower, phi.upper);
    if(length(start.phi) == 1)
  	  phi = start.phi;
    lml.phi = lMarLik.phi.normAR(Xtilde, phi, eta, nu0, trend);
    current = list(lml.phi = lml.phi, phi = phi, sigmasq = NULL, alpha = NULL, mu = NULL, change.phi = 0);  	
  
    
  ## Start MCMC
    for(it in 1:iter){  	
      ## update parameters
  	  current = up.parameters.normAR(Xtilde, meta, eta, current, nu0, phi.lower, phi.upper, trend, sd.xi);
  	  	
      ## save Markov chains    
      if(it %% thin == 0){
  	    Phi[it / thin + 1] = current$phi;
  	    Sigmasq[it / thin] = current$sigmasq;
  	    Mu[it / thin, ] = current$mu + mu0;
  	    if(trend == TRUE)
  	      Alpha[it / thin ] = current$alpha;
  	      
        ## show: x0% completed
        if( (it * 10) %% iter == 0 && it != iter)
          cat(paste( (it * 100) / iter), '% completed...\n', sep = '');
        if( it == iter){
          cat(paste( (it * 100) / iter), '% completed.\n', sep = '');
          cat('\n');
        }
      } 
    }
  change.phi = current$change.phi / iter;
  }
  
  ## track time
  if(track.time == TRUE){
    t.finish = proc.time();
    cat('Time used (in second): \n')
    print(t.finish - t.start);
    cat('\n');
  }
  
  input.parameters = list(eta = eta, iter = iter, thin = thin, trend = trend, EB = EB, mu0 = mu0, nu0 = nu0, phi.lower = phi.lower, phi.upper = phi.upper, sd.xi = sd.xi, start.phi = start.phi, start.year = start.year);
  
  ## point estimates
  keep = round(iter / thin * burnin) : (iter / thin);
  if(EB == TRUE){
  	phi.est = Phi;
  	sigmasq.est = Sigmasq;
  }
  if(EB == FALSE){
  	phi.est = mean(Phi[keep]);
  	sigmasq.est = mean(Sigmasq[keep]);
  }
  if(trend == TRUE) 
  	alpha.est = mean(Alpha[keep]);
  mu.est = apply(Mu[keep, ], 2, mean);
  
  ## show summary
  if(show.summary == TRUE){
    cat('Estimates of parameters: \n')
  	cat('phi: ', round(phi.est, 4), '\n', sep = '');
  	cat('sigma2: ', round(sigmasq.est, 4), '\n', sep = '');
  	if(trend == TRUE) 
  	  cat('alpha: ', round(alpha.est, 4), '\n', sep = '');
  	cat('mu:', round(mu.est, 4), '\n', sep = ' ');
  }
    
  if(trend == TRUE)
    return( list(Phi = Phi, Sigmasq = Sigmasq, Alpha = Alpha, Mu = Mu, phi.est = phi.est, sigmasq.est = sigmasq.est, alpha.est = alpha.est, mu.est = mu.est, X = X, meta = meta, input.parameters = input.parameters, change.phi = change.phi) );
  if(trend == FALSE)
    return( list(Phi = Phi, Sigmasq = Sigmasq, Alpha = 0, Mu = Mu, phi.est = phi.est, sigmasq.est = sigmasq.est, alpha.est = 0, mu.est = mu.est, X = X, meta = meta, input.parameters = input.parameters, change.phi = change.phi) );

}



