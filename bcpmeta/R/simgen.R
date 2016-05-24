#######################################################
#### Li, Y. and Lund, R (2014)                     ####
#### Multiple Changepoint Detection Using Metadata ####
#### Section 4: simulation. Data Generator.        ####
#######################################################

simgen = function(scenario, N = 1000){
	
  n = 200; 		## length of time series
  phi.true = 0.2;	## true value of phi
  sigma.true = sqrt(0.025);	## true value of sigma
  X = NULL;
  meta = c(25, 50, 85, 125, 170, 185);
 
  ## true value of mu (normal means) in different scenarios
  if(scenario == 1)
    mu.true = rep(0, n); ## 
  if(scenario == 2)
    mu.true = c(rep(0, 49), rep(0.2, 50), rep(0.4, 50), rep(0.6, 51)); ## true value of mu 
  if(scenario == 3)
    mu.true = c(rep(0, 24), rep(-0.2, 50), rep(0.2, 25), rep(0, 101)); ## true value of mu (normal means)

  ## set seed
  set.seed(1);
    
  ## generate N independent series
  for(r in 1:N){
    z.true = rnorm(n, 0, sigma.true);	## true values of white noises	
    epsilon.true = rep(NA, n);			## true values of AR(1) errors
    epsilon.true[1] = z.true[1];
    for(i in 2:n){
      epsilon.true[i] = phi.true * epsilon.true[i - 1] + z.true[i];
    }
    X = rbind(X, mu.true + epsilon.true);		## reponse 
  }
  
  return(list(X = X, meta = meta, scenario = scenario, N = N));
}