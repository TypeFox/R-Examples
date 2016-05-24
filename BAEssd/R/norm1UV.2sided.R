norm1UV.2sided <-
function(theta0,prob,mu,scale,shape,rate){
  cat("\nLoading the 'norm1UV.2sided' suite...",
      "\n  This suite contains functions pertaining to one-sample experiment",
      "\n  involving a normally distributed response with unknown variance.",
      "\n  The hypothesis of interest has a two-sided alternative.\n\n")
  
  ### Create functions
  # function: logm
  # purpose: returns a list giving the log marginal distribution, the log 
  #          marginal under the null, and the log marginal under the alternative 
  #          hypothesis.
  logm <- function(xbar,s2,n,theta0,prob,mu,scale,shape,rate){
    # Error checks
    n <- as.integer(n)
    if(prob<=0 || prob>=1) stop("prob must be in (0,1).")
    if(sum(s2<=0)>0 | scale<=0 | shape<=0 | rate<=0){
      stop("s2, scale, shape and rate must be positive.")
    }
    if(max(c(length(n),length(theta0),length(prob),length(mu),length(scale),
            length(shape),length(rate)))>1){
      stop("n, theta0, prob, mu, scale, shape and rate should have length 1.")
    }
    
    # Change parameters
    shape <- 2*shape
    rate <- 2*rate
    
    # Log marginal under null
    z <- (xbar-theta0)/sqrt((s2+rate)/(n*(n+shape-1)))
    t.xbar <- dt(z,df=(n+shape-1))/sqrt((s2+rate)/(n*(n+shape-1)))
    gg.s2 <- dggamma(s2,shape/2,rate,0.5*(n-1))
    
    logm0 <- log(t.xbar) + log(gg.s2)
    
    # Log marginal under alternative
    z <- (xbar-mu)/sqrt((s2+rate)*(1/n+(scale^2))/(n+shape-1))
    t.xbar <- dt(z,df=(n+shape-1))/sqrt((s2+rate)*(1/n+(scale^2))/(n+shape-1))
    
    logm1 <- log(t.xbar) + log(gg.s2)
    
    # Log marginal
    logm <- log(prob*exp(logm0) + (1-prob)*exp(logm1))
    
    out <- list(logm0=logm0,logm1=logm1,logm=logm)
    return(out)
  }
  
  # function: logbf
  # purpose: compute the log bayes factor.
  logbf <- function(xbar,s2,n,theta0,prob,mu,scale,shape,rate){
    # Obtain log marginals
    marg <- logm(xbar=xbar,s2=s2,n=n,theta0=theta0,prob=prob,mu=mu,scale=scale,
        shape=shape,rate=rate)
    
    # BF
    out <- marg$logm1 - marg$logm0
    
    return(out)
  }
  
  # function: prior
  # purpose: returns the prior density.
  prior <- function(theta,sigma2,theta0,prob,mu,scale,shape,rate){
    # Error checks
    if(prob<=0 || prob>=1) stop("prob must be in (0,1).")
    if(scale<=0 | shape<=0 | rate<=0){
      stop("scale, shape, and rate must be positive.")
    }
    if(max(c(length(theta0),length(prob),length(mu),length(scale),
            length(shape),length(rate)))>1){
      stop("theta0, prob, mu, scale, shape, and rate should have length 1.")
    }
    
    # Change parameters
    shape <- 2*shape
    rate <- 2*rate
    
    # Prior
    out <- prob*(theta==theta0) + 
        (1-prob)*(theta!=theta0)*dnorm(theta,mu,(scale*sqrt(sigma2)))
    out <- out*dgamma(1/sigma2,shape=shape,rate=rate)
    
    return(out)
  }
  
  # function: post
  # purpose: returns the posterior density.
  post <- function(theta,sigma2,xbar,s2,n,theta0,prob,mu,scale,shape,rate){
    # Error checks
    n <- as.integer(n)
    if(prob<=0 || prob>=1) stop("prob must be in (0,1).")
    if(scale<=0 | shape<=0 | rate<=0){
      stop("scale, shape, and rate must be positive.")
    }
    if(max(c(length(theta0),length(prob),length(mu),length(scale),
            length(shape),length(rate),length(xbar),length(s2),length(n)))>1){
      stop("xbar,s2,n,theta0,prob,mu,scale,shape,rate should have length 1.")
    }
    
    # Change parameters
    shape <- 2*shape
    rate <- 2*rate
    
    # Joint Density
    joint <- prior(theta=theta,sigma2=sigma2,theta0=theta0,prob=prob,mu=mu,
            scale=scale,shape=shape,rate=rate)*
        dnorm(xbar,mean=theta,sd=sqrt(sigma2/n))*
        dgamma(s2,shape=(0.5*(n-1)),rate=(1/(2*sigma2)))
    
    # Posterior
    out <- joint/exp(logm(xbar,s2,n,theta0,prob,mu,scale,shape,rate)$logm)
    
    return(out)
  }
  
  
  # function: ssd.norm1UV.2sided
  # purpose: Sample size determination via the Bayesian average error based
  #            approach for this specific example.
  #
  # parameters:
  #   alpha      Scalar. Bound on the total error. Sample size will be chosen
  #                such that total error is not greater than alpha.
  #   w          Scalar. The weight to be given to Average Type-I Error when
  #                minimizing total weighted error. Larger values of w control
  #                Type-I error rates more.
  #  
  #   m          Scalar. Number of MC replicates to generate to perform
  #                integration (default=2500).
  #   minn       Scalar. The minimum sample size to consider (default=2).
  #   maxn       Scalar. The maximum sample size to consider (default=1000).
  #   all        Boolean. If FALSE (default), the function returns when the 
  #                minimum sample size is determined. If TRUE, all sample sizes
  #                in the range of [minn,maxn] are considered. This is useful
  #                for tracing out the total error as the sample size increases
  #
  #   additional parameters correspond to this particular set-up.
  ssd.norm1UV.2sided <- function(alpha,w,theta0,prob,mu,scale,shape,rate,
      m=2500,minn=3,maxn=1000,all=FALSE){
    
    ### Error checking
    minn <- max(floor(minn),2)
    maxn <- max(floor(maxn),2)
    if(maxn < minn) stop("Minimum n greater than maximum n.")
    if(w <= 0 | w>=1) stop("Weight w must be in (0,1).")
    if(alpha<=0) stop ("The bound alpha must be positive.")
    if(prob<=0 || prob>=1) stop("prob must be in (0,1).")
    if(scale<=0 | shape<=0 | rate<=0){
      stop("scale, shape and rate must be positive.")
    }
    if(max(c(length(theta0),length(prob),length(mu),length(scale),
            length(shape),length(rate)))>1){
      stop("theta0, prob, mu, scale, shape and rate should have length 1.")
    }
    
    # Change parameters
    shape <- 2*shape
    rate <- 2*rate
    
    ### Define necessary internal functions
    # function: AE1
    # purpose: Create Average Type-I Error with BF as Statistic
    AE1 <- function(t,n,m,theta0,prob,mu,scale,shape,rate){
      # Generate data
      gg.s2 <- rggamma(m,shape/2,rate,0.5*(n-1))
      t.xbar <- rt(m,df=(n+shape-1))*sqrt((gg.s2+rate)/(n*(n+shape-1))) + theta0
      
      out <- mean(logbf(t.xbar,gg.s2,n,theta0,prob,mu,scale,shape,rate)>t)
      return(out)
    }
    
    # function: AE2
    # purpose: Create Average Type-II Error with BF as Statistic
    AE2 <- function(t,n,m,theta0,prob,mu,scale,shape,rate){
      # Generate data
      gg.s2 <- rggamma(m,shape/2,rate,0.5*(n-1))
      t.xbar <- rt(m,df=(n+shape-1))*
          sqrt((gg.s2+rate)*(1/n+(scale^2))/(n+shape-1)) + mu
      
      out <- mean(logbf(t.xbar,gg.s2,n,theta0,prob,mu,scale,shape,rate)<=t)
      return(out)
    }
    
    
    ### Set-up
    n <- minn
    history <- data.frame(n=NA,AE1=NA,AE2=NA,TWE=NA,TE=NA)
    
    ### Iterate process
    repeat{        
      # Get Errors, check if criteria met
      err1 <- AE1(t=log(w/(1-w)),n=n,m=m,theta0=theta0,prob=prob,mu=mu,
          scale=scale,shape=shape,rate=rate)
      err2 <- AE2(t=log(w/(1-w)),n=n,m=m,theta0=theta0,prob=prob,mu=mu,
          scale=scale,shape=shape,rate=rate)
      TWE <- w*err1 + (1-w)*err2
      TE <- err1 + err2
      
      # Record results
      history[n-minn+1,] <- c(n,err1,err2,TWE,TE)
      
      # Determine if n satisfies criteria, and stop if applicable
      if(((TE<=alpha) && !all) || n==maxn) break
      
      n <- n+1
    }
    
    # Determine selected sample size
    if(sum(history$TE<=alpha)>0) n <- min(history$n[history$TE<=alpha])
    if(sum(history$TE<=alpha)==0) n <- maxn
    attributes(n) <- list(alpha=alpha,w=w,TE=history$TE[history$n==n])
    
    out <- list(call=sys.call(),history=history,n=n)
    class(out) <- "BAEssd"
    
    return(out)
  }
  
  
  ### Assign defaults  
  if(!missing(theta0)) formals(logm)$theta0 <- theta0
  if(!missing(prob)) formals(logm)$prob <- prob
  if(!missing(mu)) formals(logm)$mu <- mu
  if(!missing(scale)) formals(logm)$scale <- scale
  if(!missing(shape)) formals(logm)$shape <- shape
  if(!missing(rate)) formals(logm)$rate <- rate
  
  formals(logbf) <- formals(logm)
  formals(prior)[-c(1:2)] <- formals(logm)[-c(1:3)]
  formals(post)[-c(1:2)] <- formals(logm)
  formals(ssd.norm1UV.2sided)[c(3:8)] <- formals(logm)[-c(1:3)] 
  
  out <- list(logm=logm,logbf=logbf,prior=prior,post=post,
      ssd.norm1UV.2sided=ssd.norm1UV.2sided)
  return(out)
}
