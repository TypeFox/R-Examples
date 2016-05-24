binom1.1sided <-
function(p0,a,b){
  cat("\nLoading the 'binom1.1sided' suite...",
      "\n  This suite contains functions pertaining to a one sample experiment",
      "\n  with a binary outcome. The hypothesis of interest has a one-sided",
      "\n  alternative.\n\n")
  
  ### Create functions
  # function: logm
  # purpose: returns a list giving the log marginal distribution, the log 
  #          marginal under the null, and the log marginal under the alternative 
  #          hypothesis.
  logm <- function(x,n,p0,a,b){
    # Error checks
    x <- as.integer(x)
    n <- as.integer(n)
    if(p0<=0 || p0>=1) stop("p0 must be in  (0,1).")
    if(a<=0 || b<=0) stop("a,b must be positive.")
    if(max(c(length(n),length(p0),length(a),length(b)))>1){
      stop("n, p0, a, and b should have length 1.")
    }
    
    # Log marginal
    logm <- lchoose(n,x) + lbeta(a+x,b+n-x) - lbeta(a,b)
    
    # Prior and posterior probabilities of null
    lF.prior <- pbeta(p0,shape1=a,shape2=b,log.p=TRUE)
    lF.post <- pbeta(p0,shape1=a+x,shape2=b+n-x,log.p=TRUE)
    
    # Log marginal under null
    logm0 <- logm + lF.post - lF.prior
    
    # Log marginal under alternative
    logm1 <- logm + log(1-exp(lF.post)) - log(1-exp(lF.prior))
    
    out <- list(logm0=logm0,logm1=logm1,logm=logm)
    return(out)
  }
  
  # function: logbf
  # purpose: compute the log bayes factor.
  logbf <- function(x,n,p0,a,b){
    # Obtain log marginals
    marg <- logm(x=x,n=n,p0=p0,a=a,b=b)
    
    # BF
    out <- marg$logm1 - marg$logm0
    
    return(out)
  }
  
  # function: prior
  # purpose: returns the prior density.
  prior <- function(p,a,b){
    # Error checks
    if(a<=0 || b<=0) stop("a,b must be positive.")
    if(max(c(length(a),length(b)))>1) stop("a and b should have length 1.")
    
    # Prior
    out <- dbeta(p,shape1=a,shape2=b)
    
    return(out)
  }
  
  # function: post
  # purpose: returns the posterior density.
  post <- function(p,x,n,a,b){
    # Error checks
    x <- as.integer(x)
    n <- as.integer(n) 
    if(a<=0 || b<=0) stop("a,b must be positive.")
    if(max(c(length(x),length(n),length(a),length(b)))>1){
      stop("x, n, a, and b should have length 1.")
    }
    
    # Posterior
    out <- dbeta(p,shape1=a+x,shape2=b+n-x)
    
    return(out)
  }
  
  ### Assign defaults
  if(!missing(p0)) formals(logm)$p0 <- p0
  if(!missing(a)) formals(logm)$a <- a
  if(!missing(b)) formals(logm)$b <- b
  
  formals(logbf) <- formals(logm)
  formals(prior)[-1] <- formals(logm)[-c(1:3)]
  formals(post)[-1] <- formals(logm)[-3]
  
  out <- list(logm=logm,logbf=logbf,prior=prior,post=post)
  return(out)
}
