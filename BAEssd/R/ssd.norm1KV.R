ssd.norm1KV <-
function(alpha,w,logm,
    minn=2,maxn=1000,all=FALSE,...){
  
  ### Error checking
  minn <- max(floor(minn),2)
  maxn <- max(floor(maxn),2)
  if(maxn < minn) stop("Minimum n greater than maximum n.")
  if(w <= 0 | w>=1) stop("Weight w must be in (0,1).")
  if(alpha<=0) stop ("The bound alpha must be positive.")
  
  ### Define necessary internal functions
  # function: AE1
  # purpose: Create Average Type-I Error with BF as Statistic
  AE1 <- function(t,n){
    f <- function(x,n){
      marg <- logm(x,n,...)
      Ts <- marg$logm1 - marg$logm0
      
      return((Ts>t)*exp(marg$logm0))
    }
    
    out <- integrate(f,lower=-Inf,upper=Inf,n=n)$value
    
    return(out)
  }
  
  # function: AE2
  # purpose: Create Average Type-II Error with BF as Statistic
  AE2 <- function(t,n){
    f <- function(x,n){
      marg <- logm(x,n,...)
      Ts <- marg$logm1 - marg$logm0
      
      return((Ts<=t)*exp(marg$logm1))
    }
    
    out <- integrate(f,lower=-Inf,upper=Inf,n=n)$value
    
    return(out)
  }
  
  
  ### Set-up
  n <- minn
  history <- data.frame(n=NA,AE1=NA,AE2=NA,TWE=NA,TE=NA)
  
  ### Iterate process
  repeat{        
    # Get Errors, check if criteria met
    err1 <- AE1(t=log(w/(1-w)),n=n)
    err2 <- AE2(t=log(w/(1-w)),n=n)
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
