#==============================================================================
# Initialize the Gaussian-Student mixture model
#
# Date: 
#   Revised: February 15, 2015
#
# Author: Jiarui Ding <jiaruid@cs.ubc.ca>
#   Department of Computer Science, UBC
#   Department of Molecular Oncology, BC Cancer Agency 
#
#' @importFrom stats kmeans aggregate sd
MixGaussInit = function (x, lambda=NULL, mu=NULL, sigma=NULL, 
                         sigma.equal=FALSE, K=2, prior=NULL) {
  n = length(x)
  x = sort(x)
  
  # Assign data to clusters by kmeans
  K.cluster = kmeans(x, K)
  if (is.null(sigma)) {
    sigma = aggregate(x, by=list(K.cluster$cluster), FUN=sd)[,2]
    
    if (sigma.equal) {
      sigma = rep(mean(sigma), K)
    }
  }
  if (is.null(mu)) {
    mu = K.cluster$centers[, 1]
  }
  if (is.null(lambda)) {
    lambda = K.cluster$size / n
  } else if (length(lambda) == 1) {
    lambda = rep(lambda, length.out=K)
    lambda = lambda/sum(lambda)
  }
  
  list(lambda=lambda, mu=mu, sigma=sigma, K=K)
}

MixStudentInit = function (x, lambda=NULL, mu=NULL, sigma=NULL, 
                           sigma.equal=FALSE, K=2, nu=NULL, prior=NULL) {
  model = MixGaussInit(x, lambda, mu, sigma, sigma.equal, K, prior)
  
  if(is.null(nu)) {
    nu = rep(5, K)
  }
  model$nu = nu
  
  return(model)
}
