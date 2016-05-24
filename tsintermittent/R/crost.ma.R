crost.ma <- function(data,h=10,w=NULL,nop=c(2,1),type=c("croston","sba","sbj"),
                     cost=c("mar","msr","mae","mse"),outplot=c(FALSE,TRUE),
                     na.rm=c(FALSE,TRUE)){
# Croston style moving averages and variants
#
# Inputs:
#   data        Intermittent demand time series.
#   h           Forecast horizon.
#   w           Moving average order. If w == NULL then moving average orders are optimised.
#               If w is a single value then the same order is used for smoothing both the 
#               demand and the intervals. If two values are provided then the second 
#               is used to smooth the intervals. 
#   nop         Specifies the number of model parameters. Used only if they are optimised.
#               1 - Demand and interval moving average order are the same
#               2 - Different demand and interval orders
#   type        Croston's method variant:
#                 1 - "croston" Croston's method;
#                 2 - "sba" Syntetos-Boylan approximation;
#                 3 - "sbj" Shale-Boylan-Johnston.
#   cost        Cost function used for optimisation
#                 "mar" - Mean absolute rate
#                 "msr" - Mean squared rate
#                 "mae" - Mean absolute error
#                 "mse" - Mean squared error
#   outplot     If TRUE a plot of the forecast is provided.
#   na.rm       A logical value indicating whether NA values should be remove using the method.
#
# Outputs:
#   model       Type of model fitted.
#   frc.in      In-sample demand rate. 
#   frc.out     Out-of-sample demand rate.
#   order       Moving averages orders for demand and interval.
#   component   List of c.in and c.out containing the non-zero demand and interval vectors for 
#               in- and out-of-sample respectively. Third element is the coefficient used to scale
#               demand rate for sba and sbj.
#
# Example:
#   crost.ma(ts.data1,outplot=TRUE)
# 
# Notes:
# Optimisation cost functions and properties described in:
# N. Kourentzes, 2014, On intermittent demand model optimisation and selection, 
# International Journal of Production Economics, 156: 180-190. 
# http://dx.doi.org/10.1016/j.ijpe.2014.06.007
# http://kourentzes.com/forecasting/2014/06/11/on-intermittent-demand-model-optimisation-and-selection/
#
# Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>

  # Defaults
  type <- tolower(type[1])
  cost <- cost[1]
  nop <- nop[1]
  outplot <- outplot[1]
  na.rm <- na.rm[1]
  
  # Prepare data
  if (class(data)=="data.frame"){
    if (ncol(data)>1){
      warning("Data frame with more than one columns. Using only first one.")
    }
    data <- data[[1]]
  }
  if (na.rm == TRUE){
    data <- data[!is.na(data)]
  }
  n <- length(data)
  
  # Croston decomposition
  nzd <- which(data != 0)               # Find location on non-zero demand
  k <- length(nzd)
  z <- data[nzd]                        # Demand
  x <- c(nzd[1],nzd[2:k]-nzd[1:(k-1)])  # Intervals
  
  # Optimise parameters if requested
  if (is.null(w)){
    w <- crost.ma.opt(data,type,cost,nop,k-1)
  }
  
  # Assign parameters
  if (length(w)==1){
    k.demand <- w[1]
    k.interval <- w[1]
  } else {
    k.demand <- w[1]
    k.interval <- w[2]
  }
  
  # Set coefficient
  if(type == "sba"){
    coeff <- k.interval/(k.interval+1)
  } else if(type == "sbj"){
    coeff <- (k.interval-1)/k.interval
  } else {
    coeff <- 1
  }
  
  # Calculate MA using filter
  zfit <- filter(z,rep(1/k.demand,k.demand),sides=2)
  zfit <- c(rep(NA,k.demand-1),zfit[!is.na(zfit)])
  xfit <- filter(x,rep(1/k.interval,k.interval),sides=2)
  xfit <- c(rep(NA,k.interval-1),xfit[!is.na(xfit)])
  cc <- coeff * zfit/xfit
  
  # Calculate in-sample demand rate
  frc.in <- x.in <- z.in <- rep(NA,n)
  tv <- c(nzd+1,n)  # Time vector used to create frc.in forecasts
  for (i in 1:k){
    if (tv[i]<=n){
      frc.in[tv[i]:min(c(tv[i+1],n))] <- cc[i]
      x.in[tv[i]:min(c(tv[i+1],n))] <- xfit[i]
      z.in[tv[i]:min(c(tv[i+1],n))] <- zfit[i]
    }
  }
  
  # Forecast out-of-sample demand rate
  if (h>0){
    frc.out <- rep(cc[k],h)
    x.out <- rep(xfit[k],h)
    z.out <- rep(zfit[k],h)
  } else {
    frc.out <- x.out <- z.out <- NULL
  }
  
  # Plot
  if (outplot==TRUE){
    plot(1:n,data,type="l",xlim=c(1,(n+h)),xlab="Period",ylab="",
         xaxs="i",yaxs="i",ylim=c(0,max(data)*1.1))
    lines(which(data>0),data[data>0],type="p",pch=20)
    lines(1:n,frc.in,col="red")
    lines((n+1):(n+h),frc.out,col="red",lwd=2)
  }
  
  # Prepare demand and interval vectors for output
  c.in <- array(cbind(z.in,x.in),c(n,2),dimnames=list(NULL,c("Demand","Interval")))
  if (h>0){
    c.out <- array(cbind(z.out,x.out),c(h,2),dimnames=list(NULL,c("Demand","Interval")))
  } else {
    c.out <- NULL
  }
  c.coeff <- coeff
  comp <- list(c.in=c.in,c.out=c.out,coeff=coeff)
  
  return(list(model=paste("ma.",type,sep=""),frc.in=as.numeric(frc.in),
              frc.out=frc.out,order=c(k.demand,k.interval),components=comp))
  
}

#-------------------------------------------------
crost.ma.opt <- function(data,type=c("croston","sba","sbj"),
                         cost=c("mar","msr","mae","mse"),nop=c(2,1),k){
  
# Optimisation function for Croston and variants
  
  type <- type[1]
  cost <- cost[1]
  nop <- nop[1]
  
  if (nop==1){
    jmax=1
  } else {
    jmax=k
  }
  
  err <- array(NA,c(k,jmax))
  
  # Grid search to avoid integer optimisation
  for (i in 1:k){
    for (j in 1:jmax){
      if (nop==2){
        err[i,j] <- crost.ma.cost(c(i,j),data,cost,type)
      } else {
        err[i,j] <- crost.ma.cost(c(i,i),data,cost,type)
      }
    }
  }
  
  wopt <- as.numeric(which(err==min(err,na.rm=TRUE),arr.ind=TRUE))
  
  if (nop==1){
    wopt <- wopt[1]
  }
    
  return(wopt)
  
}

#-------------------------------------------------
crost.ma.cost <- function(w,data,cost,type){
  # Cost functions for Croston and variants 
  frc.in <- crost.ma(data=data,w=w,h=0,type=type)$frc.in
  
  if (cost == "mse"){
    E <- data - frc.in  
    E <- E[!is.na(E)]
    E <- mean(E^2)
  } else if(cost == "mae"){
    E <- data - frc.in  
    E <- E[!is.na(E)]
    E <- mean(abs(E))
  } else if(cost == "mar"){
    n <- length(data)
    temp <- cumsum(data)/(1:n)
    n <- ceiling(0.3*n)
    temp[1:n] <- temp[n]
    E <- abs(frc.in - temp)
    E <- E[!is.na(E)]
    E <- sum(E)
  } else if(cost == "msr"){
    n <- length(data)
    temp <- cumsum(data)/(1:n)
    n <- ceiling(0.3*n)
    temp[1:n] <- temp[n]
    E <- (frc.in - temp)^2
    E <- E[!is.na(E)]
    E <- sum(E)    
  }

  return(E)
  
}