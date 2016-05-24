norm1KV.2sided <-
function(sigma,theta0,prob,mu,tau){
  cat("\nLoading the 'norm1KV.2sided' suite...",
      "\n  This suite contains functions pertaining to one-sample experiment",
      "\n  involving a normally distributed response with known variance. The",
      "\n  hypothesis of interest has a two-sided alternative.\n\n")
  
  ### Create functions
  # function: logm
  # purpose: returns a list giving the log marginal distribution, the log 
  #          marginal under the null, and the log marginal under the alternative 
  #          hypothesis.
  logm <- function(xbar,n,sigma,theta0,prob,mu,tau){
    # Error checks
    n <- as.integer(n)
    if(prob<=0 || prob>=1) stop("prob must be in (0,1).")
    if(sigma<=0 || tau<=0) stop("sigma and tau must be positive.")
    if(max(c(length(n),length(sigma),length(theta0),length(prob),
            length(mu),length(tau)))>1){
      stop("n, sigma, theta0, prob, mu, and tau should have length 1.")
    }
    
    # Log marginal under null
    logm0 <- dnorm(xbar,mean=theta0,sd=(sigma/sqrt(n)),log=TRUE)
    
    # Log marginal under alternative
    logm1 <- dnorm(xbar,mean=mu,sd=sqrt(sigma^2/n + tau^2),log=TRUE)
    
    # Log marginal
    logm <- log(prob*exp(logm0) + (1-prob)*exp(logm1))
    
    out <- list(logm0=logm0,logm1=logm1,logm=logm)
    return(out)
  }
  
  # function: logbf
  # purpose: compute the log bayes factor.
  logbf <- function(xbar,n,sigma,theta0,prob,mu,tau){
    # Obtain log marginals
    marg <- logm(xbar=xbar,n=n,sigma=sigma,theta0=theta0,prob=prob,mu=mu,
        tau=tau)
    
    # BF
    out <- marg$logm1 - marg$logm0
    
    return(out)
  }
  
  # function: prior
  # purpose: returns the prior density.
  prior <- function(theta,theta0,prob,mu,tau){
    # Error checks
    if(prob<=0 || prob>=1) stop("prob must be in (0,1).")
    if(tau<=0) stop("tau must be positive.")
    if(max(c(length(theta0),length(prob),length(mu),length(tau)))>1){
      stop("theta0, prob, mu, and tau should have length 1.")
    }
    
    # Prior
    out <- prob*(theta==theta0) + (1-prob)*(theta!=theta0)*dnorm(theta,mu,tau)
    
    return(out)
  }
  
  # function: post
  # purpose: returns the posterior density.
  post <- function(theta,xbar,n,sigma,theta0,prob,mu,tau){
    # Error checks
    n <- as.integer(n)
    if(prob<=0 || prob>=1) stop("prob must be in (0,1).")
    if(sigma<=0 || tau<=0) stop("sigma and tau must be positive.")
    if(max(c(length(xbar),length(n),length(sigma),length(theta0),length(prob),
            length(mu),length(tau)))>1){
      stop("xbar, n, sigma, theta0, prob, mu, and tau should have length 1.")
    }
    
    # Joint density
    joint <- prior(theta=theta,theta0=theta0,prob=prob,mu=mu,tau=tau)*
        dnorm(xbar,mean=theta,sd=(sigma/sqrt(n)))
    
    # Posterior
    out <- joint/exp(logm(xbar,n,sigma,theta0,prob,mu,tau)$logm)
    
    return(out)
  }
  
  
  # function: ssd.norm1KV.2sided
  # purpose: Fast version of ssd.norm1KV for this specific example. 
  #
  # parameters:
  #   For a description of the parameters, see ssd.norm1KV() and 
  #     norm1KV.2sided().
  ssd.norm1KV.2sided <- function(alpha,w,sigma,theta0,prob,mu,tau,
      minn=2,maxn=1000,all=FALSE){
    
    ### Error checking
    minn <- max(floor(minn),2)
    maxn <- max(floor(maxn),2)
    if(maxn < minn) stop("Minimum n greater than maximum n.")
    if(w <= 0 | w>=1) stop("Weight w must be in (0,1).")
    if(alpha<=0) stop ("The bound alpha must be positive.")
    if(prob<=0 || prob>=1) stop("prob must be in (0,1).")
    if(sigma<=0 || tau<=0) stop("sigma and tau must be positive.")
    if(max(c(length(sigma),length(theta0),length(prob),length(mu),
            length(tau)))>1){
      stop("sigma, theta0, prob, mu, and tau should have length 1.")
    }
    
    ### Define necessary internal functions
    # function: AE1
    # purpose: Create Average Type-I Error with BF as Statistic
    AE1 <- function(t,n,sigma,theta0,prob,mu,tau){
      delta <- (theta0-mu)/tau
      Q <- (2*t+delta^2-log((n^(-1)*sigma^2)/(n^(-1)*sigma^2+tau^2)))*
          (n^(-1)*sigma^2+tau^2)/(tau^2)
      
      if(Q<0) out <- 1
      if(Q>=0){      
        out <- 1 - pnorm(sqrt(Q)+((sigma*delta)/(tau*sqrt(n))),mean=0,sd=1) +
            pnorm(-sqrt(Q)+((sigma*delta)/(tau*sqrt(n))),mean=0,sd=1)
      }
      
      return(out)
    }
    
    # function: AE2
    # purpose: Create Average Type-II Error with BF as Statistic
    AE2 <- function(t,n,sigma,theta0,prob,mu,tau){
      delta <- (theta0-mu)/tau
      Q <- (2*t+delta^2-log((n^(-1)*sigma^2)/(n^(-1)*sigma^2+tau^2)))*
          (n^(-1)*sigma^2+tau^2)/(tau^2)
      
      v1 <- (sqrt(n^(-1))*sigma)/sqrt(n^(-1)*sigma^2+tau^2)
      v2 <- (tau*delta)/sqrt(n^(-1)*sigma^2+tau^2)
      
      if(Q<0) out <- 0
      if(Q>=0){
        out <- pnorm((sqrt(Q)+((sigma*delta)/(tau*sqrt(n))))*v1+v2,mean=0,sd=1)-
            pnorm((-sqrt(Q)+((sigma*delta)/(tau*sqrt(n))))*v1+v2,mean=0,sd=1)
      }
      
      return(out)
    }
    
    
    ### Set-up
    n <- minn
    history <- data.frame(n=NA,AE1=NA,AE2=NA,TWE=NA,TE=NA)
    
    ### Iterate process
    repeat{        
      # Get Errors, check if criteria met
      err1 <- AE1(t=log(w/(1-w)),n=n,sigma=sigma,theta0=theta0,prob=prob,
          mu=mu,tau=tau)
      err2 <- AE2(t=log(w/(1-w)),n=n,sigma=sigma,theta0=theta0,prob=prob,
          mu=mu,tau=tau)
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
  if(!missing(sigma)) formals(logm)$sigma <- sigma
  if(!missing(theta0)) formals(logm)$theta0 <- theta0
  if(!missing(prob)) formals(logm)$prob <- prob
  if(!missing(mu)) formals(logm)$mu <- mu
  if(!missing(tau)) formals(logm)$tau <- tau
  
  formals(logbf) <- formals(logm)
  formals(prior)[-1] <- formals(logm)[-c(1:3)]
  formals(post)[-1] <- formals(logm)
  formals(ssd.norm1KV.2sided)[c(3:7)] <- formals(logm)[-c(1:2)] 
  
  out <- list(logm=logm,logbf=logbf,prior=prior,post=post,
      ssd.norm1KV.2sided=ssd.norm1KV.2sided)
  return(out)
}
