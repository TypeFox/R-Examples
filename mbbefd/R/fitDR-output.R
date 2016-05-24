fitDR.addcomp <- function(x, theta, hessian=NULL, dist, method, convergence=0)
{
  f1 <- list(estimate=theta)
  #computes Hessian of -log Lik at estimate values
  if(all(!is.na(theta)) && method == "mle" && !is.null(hessian))
  {
    if(all(!is.na(hessian)) && qr(hessian)$rank == NCOL(hessian)){
      f1$vcov <- solve(hessian)
      f1$sd <- sqrt(diag(f1$vcov))
      f1$cor <- cov2cor(f1$vcov)
    }else
      f1$vcov <- f1$sd <- f1$cor <- NA
  }else{
    f1$vcov <- f1$sd <- f1$cor <- NA
  }
  #other fitdist components
  f1$convergence <- convergence
  f1$method <- method
  f1$n <- length(x)
  f1$data <- x
  f1$distname <- dist
  f1$weights <- f1$dots <- NULL
  f1$discrete <- FALSE
  #gof stat
  f1$loglik <- LLfunc(obs=x, theta=theta, dist=dist)
  npar <- length(theta)
  f1$aic <- -2*f1$loglik+2*npar
  f1$bic <- -2*f1$loglik+log(f1$n)*npar
  
  f1
}