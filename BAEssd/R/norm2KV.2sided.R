norm2KV.2sided <-
function(sigma,prob,mu0,tau0,mu1,tau1,mu2,tau2){
  cat("\nLoading the 'norm2KV.2sided' suite...",
      "\n  This suite contains functions pertaining to an experiment with two",
      "\n  independent samples, each with a response that is normally",
      "\n  distributed with a common known variance. The hypothesis of",
      "\n  of interest has a two-sided alternative.\n\n")
  
  require(mvtnorm)
  
  ### Create functions
  # function: logm
  # purpose: returns a list giving the log marginal distribution, the log 
  #          marginal under the null, and the log marginal under the alternative 
  #          hypothesis.
  logm <- function(xbar,n,sigma,prob,mu0,tau0,mu1,tau1,mu2,tau2){
    # Error checks
    n <- as.integer(n)
    if(sigma<=0 | tau0<=0 | tau1<=0 | tau2<=0){
      stop("sigma and taus must be positive.")
    }
    if(max(c(length(n),length(sigma),length(prob),length(mu0),length(tau0),
            length(mu1),length(tau1),length(mu2),length(tau2)))>1){
      stop("n, sigma, prob, mus, and taus should have length 1.")
    }
    
    # Restructure X
    xbar <- matrix(c(xbar),ncol=2)
    
    # Log marginal under null
    rho <- tau0^2/(sigma^2/n + tau0^2)
    Sigma <- (sigma^2/n + tau0^2)*matrix(c(1,rho,rho,1),nrow=2,ncol=2)
    logm0 <- apply(xbar,1,dmvnorm,mean=c(mu0,mu0),sigma=Sigma,log=TRUE)
    
    # Log marginal under alternative
    logm1 <- dnorm(xbar[,1],mean=mu1,sd=sqrt(sigma^2/n + tau1^2),log=TRUE) +
        dnorm(xbar[,2],mean=mu2,sd=sqrt(sigma^2/n + tau2^2),log=TRUE)
    
    # Log marginal
    logm <- log(prob*exp(logm0) + (1-prob)*exp(logm1))
    
    out <- list(logm0=logm0,logm1=logm1,logm=logm)
    return(out)
  }
  
  # function: logbf
  # purpose: compute the log bayes factor.
  logbf <- function(xbar,n,sigma,prob,mu0,tau0,mu1,tau1,mu2,tau2){
    # Obtain log marginals
    marg <- logm(xbar=xbar,n=n,sigma=sigma,prob=prob,mu0=mu0,tau0=tau0,
        mu1=mu1,tau1=tau1,mu2=mu2,tau2=tau2)
    
    # BF
    out <- marg$logm1 - marg$logm0
    
    return(out)
  }
  
  # function: prior
  # purpose: returns the prior density.
  prior <- function(theta,prob,mu0,tau0,mu1,tau1,mu2,tau2){
    # Error checks
    if(tau0<=0 | tau1<=0 | tau2<=0) stop("taus must be positive.")
    if(max(c(length(prob),length(mu0),length(tau0),length(mu1),length(tau1),
            length(mu2),length(tau2)))>1){
      stop("prob, mus, and taus should have length 1.")
    }
    
    # Restructure theta
    theta <- matrix(c(theta),ncol=2)
    
    # Prior
    out <- prob*(theta[,1]==theta[,2])*dnorm(theta[,1],mean=mu0,sd=tau0) +
        (1-prob)*(theta[,1]!=theta[,2])*dnorm(theta[,1],mean=mu1,sd=tau1)*
        dnorm(theta[,2],mean=mu2,sd=tau2)
    
    return(out)
  }
  
  # function: post
  # purpose: returns the posterior density.
  post <- function(theta,xbar,n,sigma,prob,mu0,tau0,mu1,tau1,mu2,tau2){
    # Error checks
    n <- as.integer(n)
    if(sigma<=0 | tau0<=0 | tau1<=0 | tau2<=0){
      stop("sigma and taus must be positive.")
    }
    if(max(c(length(n),length(sigma),length(prob),length(mu0),length(tau0),
            length(mu1),length(tau1),length(mu2),length(tau2)))>1){
      stop("n, sigma, prob, mus, and taus should have length 1.")
    }
    if(length(xbar)!=2) stop("xbar should have length 2.")
    
    # Restructure theta
    theta <- matrix(c(theta),ncol=2)
    xbar <- c(xbar)
    
    # Joint Density
    joint <- prior(theta=theta,prob=prob,mu0=mu0,tau0=tau0,mu1=mu1,tau1=tau1,
            mu2=mu2,tau2=tau2)*
        dnorm(xbar[1],mean=theta[,1],sd=(sigma/sqrt(n)))*
        dnorm(xbar[2],mean=theta[,2],sd=(sigma/sqrt(n)))
    
    # Posterior
    out <- joint/exp(logm(xbar,n,sigma,prob,mu0,tau0,mu1,tau1,mu2,tau2)$logm)
    
    return(out)
  }
  
  # function: ssd.norm2KV.2sided
  # purpose: Sample size determination via the Bayesian average error based
  #            approach in this particular example.
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
  ssd.norm2KV.2sided <- function(alpha,w,sigma,prob,mu0,tau0,mu1,tau1,mu2,tau2,
      m=2500,minn=2,maxn=1000,all=FALSE){
    
    ### Error checking
    minn <- max(floor(minn),2)
    maxn <- max(floor(maxn),2)
    if(maxn < minn) stop("Minimum n greater than maximum n.")
    if(w <= 0 | w>=1) stop("Weight w must be in (0,1).")
    if(alpha<=0) stop ("The bound alpha must be positive.")
    if(prob<=0 || prob>=1) stop("prob must be in (0,1).")
    if(sigma<=0 | tau0<=0 | tau1<=0 | tau2<=0){
      stop("sigma and taus must be positive.")
    }
    if(max(c(length(sigma),length(prob),length(mu0),length(tau0),
            length(mu1),length(tau1),length(mu2),length(tau2)))>1){
      stop("sigma, prob, mus, and taus should have length 1.")
    }
    
    ### Define necessary internal functions
    # function: AE1
    # purpose: Create Average Type-I Error with BF as Statistic
    AE1 <- function(t,n,m,sigma,prob,mu0,tau0,mu1,tau1,mu2,tau2){
      # Generate data
      rho <- tau0^2/(sigma^2/n + tau0^2)
      Sigma <- (sigma^2/n + tau0^2)*matrix(c(1,rho,rho,1),nrow=2,ncol=2)
      X <- rmvnorm(m,mean=c(mu0,mu0),sigma=Sigma)
      
      out <- mean(logbf(X,n,sigma,prob,mu0,tau0,mu1,tau1,mu2,tau2)>t)
      return(out)
    }
    
    # function: AE2
    # purpose: Create Average Type-II Error with BF as Statistic
    AE2 <- function(t,n,m,sigma,prob,mu0,tau0,mu1,tau1,mu2,tau2){
      # Generate data
      X <- rnorm(m,mean=mu1,sd=sqrt(sigma^2/n + tau1^2))*
          rnorm(m,mean=mu2,sd=sqrt(sigma^2/n + tau2^2))
      
      out <- mean(logbf(X,n,sigma,prob,mu0,tau0,mu1,tau1,mu2,tau2)<=t)
      return(out)
    }
    
    
    ### Set-up
    n <- minn
    history <- data.frame(n=NA,AE1=NA,AE2=NA,TWE=NA,TE=NA)
    
    ### Iterate process
    repeat{        
      # Get Errors, check if criteria met
      err1 <- AE1(t=log(w/(1-w)),n=n,m=m,sigma=sigma,prob=prob,mu0=mu0,
          tau0=tau0,mu1=mu1,tau1=tau1,mu2=mu2,tau2=tau2)
      err2 <- AE2(t=log(w/(1-w)),n=n,m=m,sigma=sigma,prob=prob,mu0=mu0,
          tau0=tau0,mu1=mu1,tau1=tau1,mu2=mu2,tau2=tau2)
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
  if(!missing(prob)) formals(logm)$prob <- prob
  if(!missing(mu0)) formals(logm)$mu0 <- mu0
  if(!missing(tau0)) formals(logm)$tau0 <- tau0
  if(!missing(mu1)) formals(logm)$mu1 <- mu1
  if(!missing(tau1)) formals(logm)$tau1 <- tau1
  if(!missing(mu2)) formals(logm)$mu2 <- mu2
  if(!missing(tau2)) formals(logm)$tau2 <- tau2
  
  formals(logbf) <- formals(logm)
  formals(prior)[-1] <- formals(logm)[-c(1:3)]
  formals(post)[-1] <- formals(logm)
  formals(ssd.norm2KV.2sided)[3:10] <- formals(logm)[-c(1:2)]
  
  out <- list(logm=logm,logbf=logbf,prior=prior,post=post,
      ssd.norm2KV.2sided=ssd.norm2KV.2sided)
  return(out)
}
