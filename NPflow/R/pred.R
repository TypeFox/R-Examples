# Get marginal likelihood for gaussian distribution with mean and
# covariance matrix being distributed from a Normal inverse Wishart
# distribution with parameters given by S

pred <- function(z, S){
  
  mu0 <- S[["mu"]]
  kappa0 <- S[["kappa"]]
  nu0 <- S[["nu"]]
  lambda0 <- S[["lambda"]]
  
  p <- length(mu0)
  S <- lambda0*(kappa0+1)/kappa0/(nu0-p+1)
  nu <- nu0-p+1
  
  out <- ((1+t(z-mu0)%*%solve(S)%*%(z-mu0)/nu)^(-(nu+p)/2)
          *exp(lgamma((nu+p)/2))/exp(lgamma((nu)/2))
          *(det(nu*p*S))^(-0.5)
  )
  return(out)
}