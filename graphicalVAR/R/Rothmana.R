Rothmana <-
function(X, Y, lambda_beta, lambda_kappa, convergence = 1e-4, gamma = 0.5, maxit.in = 100, maxit.out = 100){
  # Algorithm 2 of Rothmana, Levinaa & Ji Zhua
  
  Nvar <- ncol(X)
  n <- nrow(X)
  beta_ridge <- beta_ridge_C(X, Y, lambda_beta)
  
  # Starting values:
  beta <- matrix(0, Nvar, Nvar)  
  
  # Algorithm:
  it <- 0

  repeat{
    it <- it + 1
    kappa <- Kappa(beta, X, Y, lambda_kappa)
    beta_old <- beta
    beta <- Beta_C(kappa, beta, X, Y, lambda_beta, convergence, maxit.in) 
    
    if (sum(abs(beta - beta_old)) < (convergence * sum(abs(beta_ridge)))){
      break
    }
    
    if (it > maxit.out){
      warning("Model did NOT converge in outer loop")
      break
    }
  }
  
  ## Compute unconstrained kappa (codes from SparseTSCGM):
  ZeroIndex <- which(kappa==0, arr.ind=TRUE) ## Select the path of zeros
  WS <-  (t(Y)%*%Y - t(t(X)%*%Y) %*% beta - t(beta) %*% t(X)%*%Y + t(beta) %*% t(X)%*%X %*% beta)/(nrow(X))
  if (nrow(ZeroIndex)==0){
    out4 <- suppressWarnings(glasso(WS, rho = 0, trace = FALSE))
  } else {
    out4 <- suppressWarnings(glasso(WS, rho = 0, zero = ZeroIndex,trace = FALSE))
  }
  lik1  = determinant( out4$wi)$modulus[1]
  lik2 <- sum(diag( out4$wi%*%WS))

  pdO = sum(sum(kappa[upper.tri(kappa,diag=FALSE)] !=0))
  pdB = sum(sum(beta !=0))
  
  LLk <-  (n/2)*(lik1-lik2) 
  LLk0 <-  (n/2)*(-lik2)
  
  EBIC <-  -2*LLk + (log(n))*(pdO +pdB) + (pdO  + pdB)*4*gamma*log(2*Nvar)
  
  ### TRANSPOSE BETA!!!
  return(list(beta=t(beta), kappa=kappa, EBIC = EBIC))
}
