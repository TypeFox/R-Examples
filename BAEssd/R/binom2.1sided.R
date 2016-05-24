binom2.1sided <-
function(a1,b1,a2,b2){
  cat("\nLoading the 'binom2.1sided' suite...",
      "\n  This suite contains functions pertaining to an experiment with two",
      "\n  independent samples involving a binary outcome. The hypothesis of",
      "\n  interest has a one-sided alternative.\n\n")
  
  ### Create functions
  # function: logm
  # purpose: returns a list giving the log marginal distribution, the log 
  #          marginal under the null, and the log marginal under the alternative 
  #          hypothesis.
  logm <- function(x,n,a1,b1,a2,b2){
    # Error checks
    x <- as.integer(x)
    n <- as.integer(n)
    if(min(c(a1,b1,a2,b2))<=0) stop("shape parameters must be positive")
    if(max(c(length(n),length(a1),length(b1),length(a2),length(b2)))>1){
      stop("n and shape parameters should have length 1.")
    }
    
    # Correct format of x
    x <- matrix(c(x),ncol=2)
    x1 <- x[,1]
    x2 <- x[,2]
    
    # Create a function needed for numerical integration
    fct <- function(u,a1,b1,a2,b2) pbeta(u,a1,b1)*dbeta(u,a2,b2)
    
    # Determine prior probability of null hypothesis
    if(a1==a2 && b1==b2) prob <- 0.5
    if(a1!=a2 || b1!=b2){
      prob <- integrate(fct,lower=0,upper=1,a1=a1,b1=b1,a2=a2,b2=b2)$value
    }
    
    # Determine E[F]
    Ef <- apply(x,1,function(v) integrate(fct,lower=0,upper=1,
              a1=a1+v[1],b1=b1+n-v[1],a2=a2+v[2],b2=b2+n-v[2])$value)
    Ef <- pmax(pmin(Ef,1),0)
    
    # Log marginal under null
    logm0 <- -log(prob)+lchoose(n,x1)+lchoose(n,x2)+lbeta(a1+x1,b1+n-x1)+
        lbeta(a2+x2,b2+n-x2)+log(Ef)-lbeta(a1,b1)-lbeta(a2,b2)
    
    # Log marginal under alternative
    logm1 <- -log(1-prob)+lchoose(n,x1)+lchoose(n,x2)+lbeta(a1+x1,b1+n-x1)+
        lbeta(a2+x2,b2+n-x2)+log(1-Ef)-lbeta(a1,b1)-lbeta(a2,b2)
    
    # Log marginal
    logm <- log(prob*exp(logm0) + (1-prob)*exp(logm1))
    
    out <- list(logm0=logm0,logm1=logm1,logm=logm)
    return(out)
  }
  
  # function: logbf
  # purpose: compute the log bayes factor.
  logbf <- function(x,n,a1,b1,a2,b2){
    # Obtain log marginals
    marg <- logm(x=x,n=n,a1=a1,b1=b1,a2=a2,b2=b2)
    
    # BF
    out <- marg$logm1 - marg$logm0
    
    return(out)
  }
  
  # function: prior
  # purpose: returns the prior density.
  prior <- function(p,a1,b1,a2,b2){
    # Error checks
    if(min(c(a1,b1,a2,b2))<=0) stop("shape parameters must be positive")
    if(max(c(length(a1),length(b1),length(a2),length(b2)))>1){
      stop("shape parameters should have length 1.")
    }
    
    # Correct format of p
    p1 <- matrix(p,ncol=2)[,1]
    p2 <- matrix(p,ncol=2)[,2]
    
    # Prior
    out <- dbeta(p1,shape1=a1,shape2=b1)*dbeta(p2,shape1=a2,shape2=b2)
    
    return(out)
  }
  
  # function: post
  # purpose: returns the posterior density.
  post <- function(p,x,n,a1,b1,a2,b2){
    # Error checks
    x <- as.integer(x)
    n <- as.integer(n)
    if(min(c(a1,b1,a2,b2))<=0) stop("shape parameters must be positive")
    if(max(c(length(n),length(a1),length(b1),length(a2),length(b2)))>1){
      stop("n and shape parameters should have length 1.")
    }
    if(length(x)!=2) stop("x should have length 2.")
    
    # Correct format of p
    p1 <- matrix(p,ncol=2)[,1]
    p2 <- matrix(p,ncol=2)[,2]
    
    # Posterior
    out <- dbeta(p1,shape1=a1+x[1],shape2=b1+n-x[1])*
        dbeta(p2,shape1=a2+x[2],shape2=b2+n-x[2])
    
    return(out)
  }
  
  
  ### Assign defaults
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
