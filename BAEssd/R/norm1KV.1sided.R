norm1KV.1sided <-
function(sigma,theta0,mu,tau){
  cat("\nLoading the 'norm1KV.1sided' suite...",
      "\n  This suite contains functions pertaining to one-sample experiment",
      "\n  involving a normally distributed response with known variance. The",
      "\n  hypothesis of interest has a one-sided alternative.\n\n")
  
  ### Create functions
  # function: logm
  # purpose: returns a list giving the log marginal distribution, the log 
  #          marginal under the null, and the log marginal under the alternative 
  #          hypothesis.
  logm <- function(xbar,n,sigma,theta0,mu,tau){
    # Error checks
    n <- as.integer(n)
    if(sigma<=0 || tau<=0) stop("sigma and tau must be positive.")
    if(max(c(length(n),length(sigma),length(theta0),length(mu),length(tau)))>1){
      stop("n, sigma, theta0, mu, and tau should have length 1.")
    }
    
    # Log marginal under null
    Zprior <- pnorm((theta0-mu)/(tau))
    Zpost <- pnorm((theta0-((xbar*tau^2 + mu*sigma^2/n)/(tau^2 + sigma^2/n)))/
            sqrt((tau^2*sigma^2/n) / (tau^2 + sigma^2/n)))
    dens <- dnorm(xbar,mean=mu,sd=sqrt(sigma^2/n + tau^2),log=TRUE)
    
    logm0 <- log(Zpost) - log(Zprior) + dens
    
    # Log marginal under alternative
    logm1 <- log(1-Zpost) - log(1-Zprior) + dens
    
    # Log marginal
    logm <- dens
    
    out <- list(logm0=logm0,logm1=logm1,logm=logm)
    return(out)
  }
  
  # function: logbf
  # purpose: compute the log bayes factor.
  logbf <- function(xbar,n,sigma,theta0,mu,tau){
    # Obtain log marginals
    marg <- logm(xbar=xbar,n=n,sigma=sigma,theta0=theta0,mu=mu,tau=tau)
    
    # BF
    out <- marg$logm1 - marg$logm0
    
    return(out)
  }
  
  # function: prior
  # purpose: returns the prior density.
  prior <- function(theta,mu,tau){
    # Error checks
    if(tau<=0) stop("tau must be positive.")
    if(max(c(length(mu),length(tau)))>1){
      stop("mu, and tau should have length 1.")
    }
    
    # Prior
    out <- dnorm(theta,mean=mu,sd=tau)
    
    return(out)
  }
  
  # function: post
  # purpose: returns the posterior density.
  post <- function(theta,xbar,n,sigma,mu,tau){
    # Error checks
    n <- as.integer(n)
    if(sigma<=0 || tau<=0) stop("sigma and tau must be positive.")
    if(max(c(length(xbar),length(n),length(sigma),length(theta0),
            length(mu),length(tau)))>1){
      stop("xbar, n, sigma, theta0, mu, and tau should have length 1.")
    }
    
    # Posterior
    out <- dnorm(theta,mean=((xbar*tau^2 + mu*sigma^2/n)/(tau^2 + sigma^2/n)),
        sqrt((tau^2*sigma^2/n) / (tau^2 + sigma^2/n)))
    
    return(out)
  }
  
  
  ### Assign defaults
  if(!missing(sigma)) formals(logm)$sigma <- sigma
  if(!missing(theta0)) formals(logm)$theta0 <- theta0
  if(!missing(mu)) formals(logm)$mu <- mu
  if(!missing(tau)) formals(logm)$tau <- tau
  
  formals(logbf) <- formals(logm)
  formals(prior)[-1] <- formals(logm)[-c(1:4)]
  formals(post)[-1] <- formals(logm)[-4]
  
  out <- list(logm=logm,logbf=logbf,prior=prior,post=post)
  return(out)
}
