#==============================================================================
MixGaussSetPrior = function(x, K, mu=NULL, sigma=NULL, alpha=NULL, 
                            kapp=NULL, dof=NULL) {
  # Set Dirichlet prior (alpha) for lambda, a Gaussian prior on the mean 
  # (conditional on the variance p(mu|sigma) ~ N(mu_p, sigma/sqrt(kapp))
  # and an inverse gamma prior on the variance (dof)
  #
  # Date: 
  #   Revised: February 28, 2015
  #
  # Author: Jiarui Ding <jiaruid@cs.ubc.ca>
  #   Department of Computer Science, UBC
  #   Department of Molecular Oncology, BC Cancer Agency 
  #
  
  prior = list()
  
  if (is.null(mu)) {
    prior$mu = rep(0, K)
  } else {
    prior$mu = mu
  }
  
  if (is.null(sigma)) {
    sigma = (x - mean(x))^2 
    sigma = sum(sigma) / length(x) / K
    prior$sigma = rep(sigma, K)
  } else {
    prior$sigma = sigma
  }
  
  if (is.null(alpha)) {
    prior$alpha = rep(2, times=K)
  } else {
    prior$alpha = alpha
  }
  
  if (is.null(kapp)) {
    prior$kapp = 0.0
  } else {
    prior$kapp = kapp
  } 
  
  if (is.null(dof)) {
    prior$dof = 3
  } else {
    prior$dof = dof
  }
    
  return(prior)
}

#==============================================================================
MixStudentSetPrior = function(x, K, mu, sigma, alpha=NULL, kapp=NULL, 
                            dof=NULL) {
  prior = MixGaussSetPrior(x, K, mu, sigma, alpha, kapp, dof)
    
  return(prior)
}
