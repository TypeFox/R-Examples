LeSagePaceExperiment <- function(n=400, beta=c(0,  1, -1), rho=0.75, ndraw=1000, 
  burn.in=200, thinning=1, m=10,
  computeMarginalEffects=TRUE, ...) {
  
  if (length(beta) != 3) stop("Currently only implemented for 3 beta parameters")
  
  # design matrix with two standard normal variates as explanatory variables
  X <- cbind(intercept=1, x=rnorm(n), y=rnorm(n))
  
  # identity matrix I_n
  I_n <- sparseMatrix(i=1:n, j=1:n, x=1)
    
  # build spatial weight matrix W from random coordinates
  W <- kNearestNeighbors(x=rnorm(n), y=rnorm(n), k=6)
  
  # create samples from epsilon using independence of distributions (rnorm()) to avoid dense matrix I_n
  eps <- rnorm(n=n, mean=0, sd=1)
  z <- solve(qr(I_n - rho * W), X %*% beta + eps)
  y <- as.vector(z >= 0)  # binary observables, 0 or 1, FALSE or TRUE

  # MCMC estimation of spatial probit
  results <- sar_probit_mcmc(y, X, W, ndraw=ndraw, burn.in=burn.in, thinning=thinning, m=m,
   computeMarginalEffects=computeMarginalEffects, ...)
  results$sd <- apply(results$B, 2, sd)
  return(results)
}
