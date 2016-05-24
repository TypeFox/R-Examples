#==============================================================================
MixGaussFitEM = 
function(x, K=2, lambda=NULL, mu=NULL, sigma=NULL, prior=NULL, 
         sigma.equal=FALSE, conv.tol=1e-5, lambda.tol=1e-5, sigma.tol=1e-5, 
         max.iter=1000, verbose=FALSE) {
  # Maximum-likelihood/Maximum-a-posteriori estimation of the parameters of a 
  # univerate finite mixture of Gasussians with possible equality constraints 
  # on the variance
  #
  # Args:
  #   x: the input vector
  #
  # Returns:
  #   The estimated model
  #
  # Date: 
  #   Revised: February 15, 2015
  #
  # Author: Jiarui Ding <jiaruid@cs.ubc.ca>
  #   Department of Computer Science, UBC
  #   Department of Molecular Oncology, BC Cancer Agency 
  #
  
  x       = as.vector(x)
  N       = length(x)

  # Initialize model and check input parameters
  model   = MixGaussInit(x=x, K=K, lambda=lambda, mu=mu, sigma=sigma, 
                         sigma.equal=sigma.equal, prior=prior)
  
  if (missing(prior)) {
    prior = MixGaussSetPrior(x=x, K=K, mu=mu, sigma=sigma)
  }
  
  lambda  = model$lambda
  mu      = model$mu
  sigma   = model$sigma  
  iter    = 1

  loglik.hist = rep(0, times=max.iter)
  expec = .Call("sexp_mix_gauss_prob", 
             s_x      = as.double(x), 
             s_lambda = as.double(lambda), 
             s_mu     = as.double(mu),
             s_sigma  = as.double(sigma))
  post.prob = expec$post
  res = expec$res
    
  prior.like     = MixGaussPriorLoglik(lambda, mu, sigma, prior, sigma.equal)
  loglik.iter    = expec$loglik + prior.like
  loglik.hist[iter] = loglik.iter
        
  done = FALSE
  while (!done) {
    r     = colSums(post.prob)
    Xbar  = colSums(post.prob*x) / r
    Sbar  = colSums(post.prob * res)
    Sbar0 = Xbar - prior$mu
    Sbar0 = Sbar0 * Sbar0
                
    lambda = r + prior$alpha - 1;
    lambda = lambda / (N + sum(prior$alpha) - K)
        
    for (k in 1:K) {
      mu[k] = r[k]*Xbar[k] + prior$kapp * prior$mu[k]
      mu[k] = mu[k] / (r[k] + prior$kapp)
    } 
      
    for (k in 1:K) {
      if (sigma.equal) {
        den  = prior$dof + N + K + 2        
        num  = prior$sigma[1]^2 + sum(Sbar + prior$kapp*r / 
          (prior$kapp+r) * Sbar0)
        
        sigma[1:K] = sqrt(num / den)
        break
      } else {
        den  = prior$dof + r[k] + 1 + 2
        num  = prior$sigma[k]^2 + Sbar[k] + prior$kapp*r[k] / 
          (prior$kapp+r[k]) * Sbar0[k]
        
        sigma[k] = sqrt( num / den )
      }
    }
    
    if (any(sigma < sigma.tol)) {
      cat(paste("WARNING! Sigma < ", sigma.tol, sep=""), "\n")
    }
    
    expec = .Call("sexp_mix_gauss_prob", 
               s_x      = as.double(x), 
               s_lambda = as.double(lambda), 
               s_mu     = as.double(mu),
               s_sigma  = as.double(sigma))
    post.prob = expec$post
    res = expec$res
    loglik.iter = expec$loglik
      
    prior.like  =  MixGaussPriorLoglik(lambda, mu, sigma, prior, sigma.equal)
    loglik.iter = loglik.iter + prior.like
    iter        = iter + 1
    loglik.hist[iter] = loglik.iter
    done = MixGaussCheckConverge(lambda, mu, sigma, loglik.hist, iter, 
                                 max.iter, conv.tol, lambda.tol, verbose)
  } 
  loglik.hist = loglik.hist[1:iter] 
  
  if (sigma.equal == TRUE) {
    num.para = (K-1) + K + 1
  } else {
    num.para = (K-1) + K + K
  }
  
  AIC = 2 * loglik.iter - 2 * num.para
  BIC = 2 * loglik.iter - num.para * log(N)
  
  colnames(post.prob) = c(paste("comp", ".", 1:K, sep = ""))
  ret = list(x=x, lambda=lambda, mu=mu, sigma=sigma, posterior=post.prob, 
             loglik=loglik.hist, AIC=AIC, BIC=BIC, type="gauss")
  return(ret)
}

