###################################################################
######    Bayesian change point detection using metadata    #######
######    MCMC in the model space: M-H and random swapping  #######
###################################################################

######## Inputs #######
## X: observations, assumed to be normal with mean shifts
## meta: 1 if meta point, 0 otherwse. same length as X. meta[1] = 0
## trend: TRUE or FALSE, whether to include linear trend with time
## parameters: mu0, nu0, arho, brho, phi.lower, phi.upper

######## Outputs ########
## map200: 200 * (n + 3) matrix
##         200 top models with highest posterior probability that MCMC visited
##         along with lpost, lml, lerr
## Eta: (iter / thin + 1) * n matrix
##         MCMC samples of change point configuration (model) eta 



####### NOTE #######
## Default strategy: use mean(X) as mu0 if not speficied 

bcpmeta.model = function(X, meta, iter = 1e4, thin = 1e1, trend = TRUE, 
					EB = TRUE, mu0 = NULL, nu0 = 5, a1 = 1, a2 = 1, 
					b1 = 19, b2 = 3, phi.lower = -0.99, phi.upper = 0.99, 
					start.eta = NULL, track.time = TRUE, show.summary = TRUE, 
					start.year = 1, meta.year = FALSE){

  t.start = proc.time();

  ## pre-processing data
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
  
  ## output matrices
  map200 = matrix(NA, ncol = n + 3, nrow = 200);
  colnames(map200) = c(paste('T', 1:n, sep = ''), 'lml', 'phi.eb', 'lpost');
  Eta = matrix(NA, ncol = n, nrow = round(iter / thin) + 1);
  colnames(Eta) = paste('T', 1:n, sep = '');

  ## initial values of eta 
  if(length(start.eta) == 0)
    eta = rbinom(n, 1, sort(runif(2, 0, 0.2))[meta + 1]);	## vector of change point, start value
  if(length(start.eta) == n)
    eta = start.eta;
  if(length(start.eta) < n && length(start.eta) > 0)
    eta = loc2dirac(start.eta, n);
  if(length(start.eta) > n)
    stop('Error in start.eta: length of changepoint configuration cannot exceed length of time series.')
  eta[1] = 0;
  
  m.out = sum(eta[(which(meta == 0))] == 1);
  m.in = sum(eta[(which(meta == 1))] == 1);
  if(EB == FALSE)
    tmp = lMarLik.normAR(Xtilde, eta, nu0, phi.lower, phi.upper, trend); 
  if(EB == TRUE)
    tmp = lMarLik.EB.normAR(Xtilde, eta, nu0, phi.lower, phi.upper, trend); 
  lpost = tmp$lml + lbeta(a1 + m.out, b1 + n - 1 - sum(meta) - m.out) + lbeta(a2 + m.in, b2 + sum(meta) - m.in) - lbeta(a1, b1) - lbeta(a2, b2);
  if(EB == FALSE)
    current = list(eta = eta, lml = tmp$lml, lerr = tmp$lerr, lpost = lpost, change.eta = TRUE);
  if(EB == TRUE)
    current = list(eta = eta, lml = tmp$lml, phi = tmp$phi.eb, lpost = lpost, change.eta = TRUE);
      
  Eta[1, ] = current$eta;  
  map200[1, ] = unlist(current)[1 : (n + 3)];

  ## Start MCMC
  for(it in 1:iter){
  	
  	## Metropolis-Hastings update gamma  	
  	current = up.eta.MH.normAR(Xtilde, meta, current, a1, a2, b1, b2, nu0, phi.lower, phi.upper, trend, EB); 
   	
  	if(current$change.eta == TRUE){
      tmpk = which(current$lpost < map200[, n + 3]);
      if(length(tmpk) == 0)
        k = 0;
      if(length(tmpk) > 0)
  	    k = max(tmpk);
  	  if(k < 200 && is.na(map200[k + 1, n + 3]))
  	    map200[(k + 1), ] = unlist(current)[1 : (n + 3)];
  	  if(k < 200 && sum(map200[k + 1, 1:n] != current$eta) > 0){
  	    if(k == 199)
  	      map200[200, ] = unlist(current)[1 : (n + 3)];
  	    if(k < 199){
  	      map200[(k + 2):200, ] = map200[(k + 1):199, ];
  	      map200[(k + 1), ] = unlist(current)[1 : (n + 3)];
  	    }
  	  }
  	}
  	
  	## simple random swapping update gamma  	
  	current = up.eta.RS.normAR(Xtilde, meta, current, a1, a2, b1, b2, nu0, phi.lower, phi.upper, trend, EB); 
  	
  	if(current$change.eta == TRUE){
      tmpk = which(current$lpost < map200[, n + 3]);
      if(length(tmpk) == 0)
        k = 0;
      if(length(tmpk) > 0)
  	    k = max(tmpk);
  	  if(k < 200 && is.na(map200[k + 1, n + 3]))
  	    map200[(k + 1), ] = unlist(current)[1 : (n + 3)];
  	  if(k < 200 && sum(map200[k + 1, 1:n] != current$eta) > 0){
  	    if(k == 199)
  	      map200[200, ] = unlist(current)[1 : (n + 3)];
  	    if(k < 199){
  	      map200[(k + 2):200, ] = map200[(k + 1):199, ];
  	      map200[(k + 1), ] = unlist(current)[1 : (n + 3)];
  	    }
  	  }
  	}
  	
    ## save Markov chains    
    if(it %% thin == 0){
  	  Eta[it / thin + 1, ] = current$eta;
      
      ## show: x0% completed
      if( (it * 10) %% iter == 0 && it != iter)
        cat(paste( (it * 100) / iter), '% completed...\n', sep = '');
      if( it == iter){
        cat(paste( (it * 100) / iter), '% completed.\n', sep = '');
        cat('\n');
      } 
    }
  }
  
  if(EB == FALSE){
  	map200[, n + 2] = NA;
  }
  
  ## track time
  if(track.time == TRUE){
    t.finish = proc.time();
    cat('Time used (in second): \n')
    print(t.finish - t.start);
    cat('\n');
  }
  
  input.parameters = list(iter = iter, thin = thin, trend = trend, EB = EB, mu0 = mu0, nu0 = nu0, a1 = a1, a2 = a2, b1 = b1, b2 = b2, phi.lower = phi.lower, phi.upper = phi.upper, start.year = start.year);
  
  ## show summary: top 5 models
  if(show.summary == TRUE){
    cat('Top 5 changepoint configurations: \n');
    for(r in 1:5)
  	  cat(which(map200[r, 1:n] == 1) + start.year - 1, '\n', sep = ' ');
  }
  
  return( list(Eta = Eta, map200 = map200, X = X, meta = meta, input.parameters = input.parameters) );
}



