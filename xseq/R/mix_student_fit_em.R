#==============================================================================
#' @importFrom stats optimize 
#' 
MixStudentFitEM = 
function(x, K=2, lambda=NULL, mu=NULL, sigma=NULL, prior=NULL, 
         sigma.equal=FALSE, conv.tol=1e-5, lambda.tol=1e-5, sigma.tol=1e-5, 
         max.iter=1000, verbose=FALSE, nu.max=50, nu.equal=FALSE) {
  # Maximum-likelihood/Maximum-a-posteriori estimation of the parameters of a 
  # univerate finite mixture of students with possible equality constraints 
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
  model   = MixStudentInit(x=x, K=K, lambda=lambda, mu=mu, sigma=sigma, 
                         sigma.equal=sigma.equal, prior=prior)
  
  if(missing(prior)) {
    prior = MixStudentSetPrior(x=x, K=K, mu=mu, sigma=sigma)
  }
  
  lambda  = model$lambda
  mu      = model$mu
  sigma   = model$sigma  
  iter    = 1
  
  nu.update = optimize(EstimateDofShare, x=x, K=K, lambda=lambda, 
                mu=mu, sigma=sigma, interval=c(0.1, nu.max))
  if (nu.equal == TRUE) {
    nu = rep(nu.update[[1]], K)
  } else {
    stop("No time to implement!")
  }
  
  # Initialization
  loglik.hist = rep(0, times=max.iter)
  expec = MixStudentProb(x, K, lambda, mu, sigma, nu)
  post.prob = expec$tau 
  
  res   = sapply(seq(1:K), function(i) (x-mu[i])^2)
  mahal = sapply(seq(1:K), function(i) res[, i] / sigma[i]^2)
  eu = sapply(seq(1:K), function(i) (nu[i] + 1) / (nu[i] + mahal[, i])) 
    
  prior.like     = MixStudentPriorLoglik(lambda, mu, sigma, 
                                         prior, sigma.equal)
  loglik.iter    = expec$loglik + prior.like
  loglik.hist[iter] = loglik.iter
        

  done = FALSE
  while (!done) {
    r     = colSums(post.prob)
    tau.eu  = sapply(seq(1:K), function(i) post.prob[,i] * eu[,i])
    tau.eu.sum = colSums(tau.eu)
                    
    lambda = r + prior$alpha - 1;
    lambda = lambda / (N + sum(prior$alpha) - K)
    
    mu = colSums(tau.eu * x) + prior$kapp * prior$mu
    mu = mu / (tau.eu.sum + prior$kapp)
    
    Xbar  = colSums(tau.eu * x) / tau.eu.sum
    Sbar  = colSums(tau.eu * res)
    
    Sbar0 = Xbar - prior$mu
    Sbar0 = Sbar0 * Sbar0
  
    for (k in 1:K) {
      if (sigma.equal) {
        den  = N + prior$dof + K + 2     
        num  = sum(Sbar + prior$kapp * tau.eu.sum / 
          (prior$kapp + tau.eu.sum) * Sbar0) + prior$sigma[1]^2
        sigma[1:K] = sqrt(num / den)
        
        break
      } else {
        stop("No time to implement!")
      }
    }
    
    nu.update = optimize(EstimateDofShare, x=x, K=K, lambda=lambda, 
                  mu=mu, sigma=sigma, interval=c(0.1, nu.max))
    if (nu.equal == TRUE) {
      nu = rep(nu.update[[1]], K)
    } else {
      stop("No time to implement!")
    }
       
    
    expec = MixStudentProb(x, K, lambda, mu, sigma, nu)
    res   = sapply(seq(1:K), function(i) (x-mu[i])^2)
    
    post.prob = expec$tau
    mahal = sapply(seq(1:K), function(i) (x-mu[i])^2 / sigma[i]^2)
    eu = sapply(seq(1:K), function(i) (nu[i] + 1) / (nu[i] + mahal[, i]))
    loglik.iter = expec$loglik    
      
    prior.like  = MixStudentPriorLoglik(lambda, mu, sigma, prior, sigma.equal)
    loglik.iter = loglik.iter + prior.like
    iter        = iter + 1
    loglik.hist[iter] = loglik.iter
    done = MixStudentCheckConverge(lambda, mu, sigma, nu, loglik.hist, 
             iter, max.iter, conv.tol, lambda.tol, verbose)
  } 
  
  if (sigma.equal == TRUE & nu.equal == TRUE) {
    num.para = (K-1) + K + 1 + 1
  } else if (sigma.equal == FALSE & nu.equal == FALSE) {
    num.para = (K-1) + K + K + K
  } else {
    num.para = (K-1) + K + K + 1
  }
  
  AIC = 2 * loglik.iter - 2 * num.para
  BIC = 2 * loglik.iter - num.para * log(N)
  
  loglik.hist = loglik.hist[1:iter] 
  colnames(post.prob) = c(paste("comp", ".", 1:K, sep = ""))
  ret = list(x=x, lambda=lambda, mu=mu, sigma=sigma, nu=nu, 
             posterior=post.prob, loglik=loglik.hist, 
             AIC=AIC, BIC=BIC, type="student")
  return(ret)
}

#==============================================================================
EstimateDofShare = function(x, K, lambda, mu, sigma, nu) {
  log.lambda = log(lambda)
  log.tau = sapply(seq(1:K), function(i) log.lambda[i] + 
                     StudentLogProb(x, mu[i], sigma[i], nu))
  nloglik = log(rowSums(exp(log.tau)))
  nloglik = -sum(nloglik)
  
  return(nloglik)
}

#==============================================================================
MixStudentProb = function(x, K, lambda, mu, sigma, nu) {  
  log.lambda = log(lambda)
  log.tau = sapply(seq(1:K), function(i) log.lambda[i] + 
                 StudentLogProb(x, mu[i], sigma[i], nu[i]))
  
  tau = exp(log.tau)
  tau.sum = rowSums(tau)
  
  loglik  = log(tau.sum)
  loglik  = sum(loglik)
  tau = tau / tau.sum
  
  return(list(tau=tau, loglik=loglik))
}

## TODO: implement in C
StudentLogProb = function(x, mu, sigma, nu) {
  x = x - mu
  sigma = sigma^2
  
  s = x^2 / (sigma)
  logp  = lgamma(nu/2 + 0.5) - lgamma(nu/2) - 0.5*log(sigma) - 
          0.5*log(nu) - 0.5*log(pi)
  
  logp = logp - (nu + 1) / 2 * log1p(s/nu)
  
  return(logp)
}




