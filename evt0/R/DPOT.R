DPOT = function (x, cov=0.01, c=0.75, th=0.1,  nd=1000)
{
  
  # Checking for plausible inputes
  if (!is.numeric(x)) 
  {
    stop("Some of the data points are not real.")
  }
  
  if ( cov < 0 || cov > 1 || !is.numeric(cov)) 
  {
    stop("coverage (cov) must lie between 0 and 1.")
  }
  
  if ( c < 0 || c > 1 || !is.numeric(c)) 
  {
    stop("q must lie between 0 and 1.")
  }
  
  if ( th < 0 ||th > 1 || !is.numeric(th)) 
  {
    stop("threshold(th) must lie between 0 and 1.")
  }
    
  if( nd < 0 || nd > length(x) || !is.numeric(nd) || nd != as.integer(nd)) 
  {
    stop("Check the number of days for computation (nd).")
  }
  
  
  
  coverage = cov
  
  #### log-likelihood function which takes three arguments: theta is
  #### the vector of parameters, y the excesses and x the durations:
  gpdlik <- function(theta,y,x){
    alpha1 <- theta[1]
    gamma <- theta[2]
    n<-length(y)
    logl<- -sum(log(alpha1*(1/x)^c))-(1/gamma+1)*sum(log(1+gamma*y/(alpha1*(1/x)^c)))
    return(-logl)
  }


  xx <- x
  a <- xx*-1 
  #a <- a[,1]
  b <- a[1:nd] #from day 1 until day nd
  len <- length(b)
  
  #### Calculation of excesses and durations
  #### since the preceding 3 excesses
  b_sort <-sort(b)
  u <-b_sort[floor(th*len)]
  bb <- b[b>u]
  bb <- bb-u
  duration <- 1
  j <- 1

  xexc <- rep(0,times=length(bb))
  
  for(ii in 1:len)
  {
    if (b[ii]>u)
    {
      xexc[j] <- duration
      duration <- 1
      j <- j+1
    }
    else 
    {
      duration <-duration+1
    }
  }
  
  lag1_xexc <-rep(0,times=length(bb))
  d2 <-rep(0,times=length(xexc))
  limit <- length(xexc)-1
  xxxx <- xexc[1:limit]
  lag1_xexc <- c(0, xxxx)
  limit2 <- length(xexc)-2
  xxxx <- xexc[1:limit2]
  lag2_xexc <- c(0, 0, xxxx)

  
  limit3 <- length(bb)
  bb <- bb[3:limit3]
  xexc <- xexc[3:limit3]
  lag1_xexc <- lag1_xexc[3:limit3]
  lag2_xexc <- lag2_xexc[3:limit3]
  d3 <- xexc+lag1_xexc+lag2_xexc
  
  #### durations since the preceding 3 excesses (v=3)
  #### We use the optim with Nelder and Mead algorithm to
  #### maximize the log likelihood
  model <- optim(c(0.5,0.5), gpdlik, y=bb, x=d3)
  mle1 <- model$par[1]
  mle2 <- model$par[2]
  
  #### With the VaR DPOT estimator we compute the forecast
  delta <- mle1*(1/(duration+xexc[length(xexc)]+xexc[length(xexc)-1]))^c
  var_forecast <- u + ((th/coverage)^mle2-1)*(delta/mle2)

  #### One-day-ahead VaR forecast:
  return(list(VaR.forecast=var_forecast, alpha = mle1, gamma = mle2))
  
}


