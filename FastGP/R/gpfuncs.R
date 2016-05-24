#R functions for Gaussian processes

ess <- function(log.lik,Y, Sig, N_mcmc,burn_in,N,flag){
  print("Running elliptical slice sampling...")
  mcmc_samples <- matrix(rep(0,(burn_in+N_mcmc)*N),nrow=N_mcmc+burn_in,byrow=T)
  if(flag==TRUE){
    norm_samples <- mvrnorm(n = burn_in+N_mcmc, mu = rep(0,N), Sigma = Sig , tol = 1e-6)
  }
  else{
    norm_samples <- rcpp_rmvnorm(n = burn_in+N_mcmc, S = Sig, mu = rep(0,N))
  }
  unif_samples <- runif(n=burn_in+N_mcmc)
  theta <- runif(n=N_mcmc+burn_in,min=0,max=2*pi)
  theta_min <- theta-2*pi
  theta_max <- theta+2*pi
  for(i in 2:(N_mcmc+burn_in)){
    f <- mcmc_samples[i-1,]
    llh_thresh <- log.lik(f,Y) + log(unif_samples[i])
    f_star <- f*cos(theta[i])+norm_samples[i,]*sin(theta[i])
    while(log.lik(f_star,Y) < llh_thresh)
    {
      if (theta[i] < 0) 
      {
        theta_min[i] <- theta[i]
      }
      else
      {
        theta_max[i] <- theta[i]
      } 
      theta[i] <- runif(n=1,min=theta_min[i],max=theta_max[i])  
      f_star <- f*cos(theta[i])+norm_samples[i,]*sin(theta[i]) 		
    }
    mcmc_samples[i,] <- f_star
  }
  return(mcmc_samples[(burn_in+1):(burn_in+N_mcmc),])
}
