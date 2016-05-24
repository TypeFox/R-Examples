binom2.2sided <-
function(prob,a0,b0,a1,b1,a2,b2){
  cat("\nLoading the 'binom2.2sided' suite...",
      "\n  This suite contains functions pertaining to an experiment with two",
      "\n  independent samples involving a binary outcome. The hypothesis of",
      "\n  interest has a two-sided alternative.\n\n")
  
  ### Create functions
  # function: logm
  # purpose: returns a list giving the log marginal distribution, the log 
  #          marginal under the null, and the log marginal under the alternative 
  #          hypothesis.
  logm <- function(x,n,prob,a0,b0,a1,b1,a2,b2){
    # Error checks
    x <- as.integer(x)
    n <- as.integer(n)
    if(prob<=0 || prob>=1) stop("prob must be in (0,1).")
    if(min(c(a0,b0,a1,b1,a2,b2))<=0) stop("shape parameters must be positive")
    if(max(c(length(n),length(prob),length(a0),length(b0),
            length(a1),length(b1),length(a2),length(b2)))>1){
      stop("n, prob, and shape parameters should have length 1.")
    }
    
    # Correct format of x
    x1 <- matrix(x,ncol=2)[,1]
    x2 <- matrix(x,ncol=2)[,2]
    
    # Log marginal under null
    logm0 <- lchoose(n,x1)+lchoose(n,x2)+lbeta(a0+x1+x2,b0+n+n-x1-x2)-
        lbeta(a0,b0)
    
    # Log marginal under alternative
    logm1 <- lchoose(n,x1)+lchoose(n,x2)+lbeta(a1+x1,b1+n-x1)+
        lbeta(a2+x2,b2+n-x2)-lbeta(a1,b1)-lbeta(a2,b2)
    
    # Log marginal
    logm <- log(prob*exp(logm0) + (1-prob)*exp(logm1))
    
    out <- list(logm0=logm0,logm1=logm1,logm=logm)
    return(out)
  }
  
  # function: logbf
  # purpose: compute the log bayes factor.
  logbf <- function(x,n,prob,a0,b0,a1,b1,a2,b2){
    # Obtain log marginals
    marg <- logm(x=x,n=n,prob=prob,a0=a0,b0=b0,a1=a1,b1=b1,a2=a2,b2=b2)
    
    # BF
    out <- marg$logm1 - marg$logm0
    
    return(out)
  }
  
  # function: prior
  # purpose: returns the prior density.
  prior <- function(p,prob,a0,b0,a1,b1,a2,b2){
    # Error checks
    if(prob<=0 || prob>=1) stop("prob must be in (0,1).")
    if(min(c(a0,b0,a1,b1,a2,b2))<=0) stop("shape parameters must be positive")
    if(max(c(length(prob),length(a0),length(b0),
            length(a1),length(b1),length(a2),length(b2)))>1){
      stop("prob and shape parameters should have length 1.")
    }
    
    # Correct format of p
    p1 <- matrix(p,ncol=2)[,1]
    p2 <- matrix(p,ncol=2)[,2]
    
    # Prior
    out <- prob*(p1==p2)*dbeta(p1,shape1=a0,shape2=b0) + (1-prob)*(p1!=p2)*
        dbeta(p1,shape1=a1,shape2=b1)*dbeta(p2,shape1=a2,shape2=b2)
    
    return(out)
  }
  
  # function: post
  # purpose: returns the posterior density.
  post <- function(p,x,n,prob,a0,b0,a1,b1,a2,b2){
    # Error checks
    x <- as.integer(x)
    n <- as.integer(n)
    if(prob<=0 || prob>=1) stop("prob must be in (0,1).")
    if(min(c(a0,b0,a1,b1,a2,b2))<=0) stop("shape parameters must be positive")
    if(max(c(length(n),length(prob),length(a0),length(b0),
            length(a1),length(b1),length(a2),length(b2)))>1){
      stop("n, prob, and shape parameters should have length 1.")
    }
    if(length(x)!=2) stop("x should have length 2.")
    
    # Joint density
    p1 <- matrix(p,ncol=2)[,1]
    p2 <- matrix(p,ncol=2)[,2]
    joint <- prior(p=p,prob=prob,a0=a0,b0=b0,a1=a1,b1=b1,a2=a2,b2=b2)*
        dbinom(x[1],size=n,prob=p1)*dbinom(x[2],size=n,prob=p2)
    
    # Posterior
    out <- joint/exp(logm(x,n,prob,a0,b0,a1,b1,a2,b2)$logm)
    
    return(out)
  }
  
  
  ### Assign defaults
  if(!missing(prob)) formals(logm)$prob <- prob
  if(!missing(a0)) formals(logm)$a0 <- a0
  if(!missing(b0)) formals(logm)$b0 <- b0
  if(!missing(a1)) formals(logm)$a1 <- a1
  if(!missing(b1)) formals(logm)$b1 <- b1
  if(!missing(a2)) formals(logm)$a2 <- a2
  if(!missing(b2)) formals(logm)$b2 <- b2
  
  formals(logbf) <- formals(logm)
  formals(prior)[-1] <- formals(logm)[-c(1:2)]
  formals(post)[-1] <- formals(logm)
  
  out <- list(logm=logm,logbf=logbf,prior=prior,post=post)
  return(out)
}
