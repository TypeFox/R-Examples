binom1.2sided <-
function(p0,prob,a,b){
  cat("\nLoading the 'binom1.2sided' suite...",
      "\n  This suite contains functions pertaining to a one sample experiment",
      "\n  with a binary outcome. The hypothesis of interest has a two-sided",
      "\n  alternative.\n\n")
  
  ### Create functions
  # function: logm
  # purpose: returns a list giving the log marginal distribution, the log 
  #          marginal under the null, and the log marginal under the alternative 
  #          hypothesis.
  logm <- function(x,n,p0,prob,a,b){
    # Error checks
    x <- as.integer(x)
    n <- as.integer(n)
    if(p0<=0 || p0>=1) stop("p0 must be in  (0,1).")
    if(prob<=0 || prob>=1) stop("prob must be in (0,1).")
    if(a<=0 || b<=0) stop("a,b must be positive.")
    if(max(c(length(n),length(p0),length(prob),length(a),length(b)))>1){
      stop("n, p0, prob, a, and b should have length 1.")
    }
    
    # Log marginal under null
    logm0 <- dbinom(x,size=n,prob=p0,log=TRUE)
    
    # Log marginal under alternative
    logm1 <- lchoose(n,x)+lbeta(a+x,b+n-x)-lbeta(a,b)
    
    # Log marginal
    logm <- log(prob*exp(logm0) + (1-prob)*exp(logm1))
    
    out <- list(logm0=logm0,logm1=logm1,logm=logm)
    return(out)
  }
  
  # function: logbf
  # purpose: compute the log bayes factor.
  logbf <- function(x,n,p0,prob,a,b){
    # Obtain log marginals
    marg <- logm(x=x,n=n,p0=p0,prob=prob,a=a,b=b)
    
    # BF
    out <- marg$logm1 - marg$logm0
    
    return(out)
  }
  
  # function: prior
  # purpose: returns the prior density.
  prior <- function(p,p0,prob,a,b){
    # Error checks
    if(p0<=0 || p0>=1) stop("p0 must be in  (0,1).")
    if(prob<=0 || prob>=1) stop("prob must be in (0,1).")
    if(a<=0 || b<=0) stop("a,b must be positive.")
    if(max(c(length(p0),length(prob),length(a),length(b)))>1){
      stop("p0, prob, a, and b should have length 1.")
    }
    
    # Prior
    out <- prob*(p==p0) + (1-prob)*(p!=p0)*dbeta(p,shape1=a,shape2=b)
    
    return(out)
  }
  
  # function: post
  # purpose: returns the posterior density.
  post <- function(p,x,n,p0,prob,a,b){
    # Error checks
    x <- as.integer(x)
    n <- as.integer(n) 
    if(p0<=0 || p0>=1) stop("p0 must be in  (0,1).")
    if(prob<=0 || prob>=1) stop("prob must be in (0,1).")
    if(a<=0 || b<=0) stop("a,b must be positive.")
    if(max(c(length(x),length(n),length(p0),length(prob),length(a),
            length(b)))>1){
      stop("x, n, p0, prob, a, and b should have length 1.")
    }
    
    # Joint density
    joint <- prior(p=p,p0=p0,prob=prob,a=a,b=b)*dbinom(x,size=n,prob=p)
    
    # Posterior
    out <- joint/exp(logm(x=x,n=n,p0=p0,prob=prob,a=a,b=b)$logm)
    
    return(out)
  }
  
  
  ### Assign defaults
  if(!missing(p0)) formals(logm)$p0 <- p0
  if(!missing(prob)) formals(logm)$prob <- prob
  if(!missing(a)) formals(logm)$a <- a
  if(!missing(b)) formals(logm)$b <- b
  
  formals(logbf) <- formals(logm)
  formals(prior)[-1] <- formals(logm)[-c(1:2)]
  formals(post)[-1] <- formals(logm)
  
  out <- list(logm=logm,logbf=logbf,prior=prior,post=post)
  return(out)
}
