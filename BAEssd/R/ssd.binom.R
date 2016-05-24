ssd.binom <-
function(alpha,w,logm,
    minn=2,maxn=1000,two.sample=FALSE,all=FALSE,...){
  
  ### Error checking
  minn <- max(floor(minn),2)
  maxn <- max(floor(maxn),2)
  if(maxn < minn) stop("Minimum n greater than maximum n.")
  if(w <= 0 | w>=1) stop("Weight w must be in (0,1).")
  if(alpha<=0) stop("The bound alpha must be positive.")
  
  ### Define necessary internal functions
  # function: AE1
  # purpose: Create Average Type-I Error
  AE1 <- function(t,Ts,logm0) sum((Ts>t)*exp(logm0))
  
  # function: AE2
  # purpose: Create Average Type-II Error
  AE2 <- function(t,Ts,logm1) sum((Ts<=t)*exp(logm1))
  
  
  ### Set-up
  n <- minn
  history <- data.frame(n=NA,AE1=NA,AE2=NA,TWE=NA,TE=NA)
  
  ### Iterate process
  repeat{
    # Create observed data
    if(!two.sample) X <- c(0:n)
    if(two.sample) X <- cbind(rep(0:n,each=(n+1)),rep(0:n,times=(n+1)))
    
    # Obtain marginal distributions
    marg <- logm(X,n,...)
    logm0 <- marg$logm0
    logm1 <- marg$logm1
    
    # Calculate logBF, the test statistic
    Ts <- logm1-logm0
    
    # Get Errors, check if criteria met
    err1 <- AE1(t=log(w/(1-w)),Ts=Ts,logm0=logm0)
    err2 <- AE2(t=log(w/(1-w)),Ts=Ts,logm1=logm1)
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
