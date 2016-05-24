crost <- function(data,h=10,w=NULL,init=c("mean","naive"),nop=c(2,1),
                  type=c("croston","sba","sbj"),cost=c("mar","msr","mae","mse"),
                  init.opt=c(TRUE,FALSE),outplot=c(FALSE,TRUE),opt.on=c(FALSE,TRUE),
                  na.rm=c(FALSE,TRUE)){
# Croston method and variants
#
# Inputs:
#   data        Intermittent demand time series.
#   h           Forecast horizon.
#   w           Smoothing parameters. If w == NULL then parameters are optimised.
#               If w is a single parameter then the same is used for smoothing both the 
#               demand and the intervals. If two parameters are provided then the second 
#               is used to smooth the intervals. 
#   init        Initial values for demand and intervals. This can be:
#                 c(z,x)  - Vector of two scalars, where first is initial demand and 
#                           second is initial interval;
#                 "naive" - Initial demand is first non-zero demand and initial interval
#                           is first interval;
#                 "mean"  - Same as "naive", but initial interval is the mean of all 
#                           in sample intervals.
#   nop         Specifies the number of model parameters. Used only if they are optimised.
#               1 - Demand and interval parameters are the same
#               2 - Different demand and interval parameters
#   type        Croston's method variant:
#                 1 - "croston" Croston's method;
#                 2 - "sba" Syntetos-Boylan approximation;
#                 3 - "sbj" Shale-Boylan-Johnston.
#   cost        Cost function used for optimisation
#                 "mar" - Mean absolute rate
#                 "msr" - Mean squared rate
#                 "mae" - Mean absolute error
#                 "mse" - Mean squared error
#   init.opt    If init.opt==TRUE then initial values are optimised. 
#   outplot     If TRUE a plot of the forecast is provided.
#   opt.on      This is meant to use only by the optimisation function. When opt.on is 
#               TRUE then no checks on inputs are performed. 
#   na.rm       A logical value indicating whether NA values should be remove using the method.
#
# Outputs:
#   model       Type of model fitted.
#   frc.in      In-sample demand rate. 
#   frc.out     Out-of-sample demand rate.
#   weights     Smoothing parameters for demand and interval.
#   initial     Initialisation values for demand and interval smoothing.
#   component   List of c.in and c.out containing the non-zero demand and interval vectors for 
#               in- and out-of-sample respectively. Third element is the coefficient used to scale
#               demand rate for sba and sbj.
#
# Example:
#   crost(ts.data1,outplot=TRUE)
#
# Notes:
# Optimisation of the methods described in:
# N. Kourentzes, 2014, On intermittent demand model optimisation and selection, 
# International Journal of Production Economics, 156: 180-190. 
# http://dx.doi.org/10.1016/j.ijpe.2014.06.007
# http://kourentzes.com/forecasting/2014/06/11/on-intermittent-demand-model-optimisation-and-selection/
#
# Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>
  
  # Defaults
  type <- tolower(type[1])
  cost <- cost[1]
  init.opt <- init.opt[1]
  outplot <- outplot[1]
  opt.on <- opt.on[1]
  na.rm <- na.rm[1]
  if (!is.numeric(init)){
    init <- init[1]
  } else {
    if (length(init>=2)){
      init <- init[1:2]
    } else {
      init <- "mean"
    }
  }
  
  # Make sure that nop is of correct lenght
  if (nop>2 || nop<1){
    nop <- 2
    warning("nop can be either 1 or 2. Overriden to 2.")
  }
  
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
  
  # Check number of non-zero values - need to have at least two
  if (sum(data!=0)<2){
    stop("Need at least two non-zero values to model time series.")
  }
  
  # Croston decomposition
  nzd <- which(data != 0)               # Find location on non-zero demand
  k <- length(nzd)
  z <- data[nzd]                        # Demand
  x <- c(nzd[1],nzd[2:k]-nzd[1:(k-1)])  # Intervals
  
  # Initialise
  if (!(is.numeric(init) && length(init)==2)){
    if (init=="mean"){
      init <- c(z[1],mean(x))
    } else {
      init <- c(z[1],x[1])
    }
  }

  # Optimise parameters if requested
  if (opt.on == FALSE){
    if (is.null(w) || init.opt == TRUE){
      wopt <- crost.opt(data,type,cost,w,nop,init,init.opt)
      w <- wopt$w
      init <- wopt$init
    }
  }
  
  # Pre-allocate memory
  zfit <- vector("numeric",k)
  xfit <- vector("numeric",k)
  
  # Assign initial values and parameters
  if (opt.on == FALSE){
    if (init[1]<0){
      stop("Initial demand cannot be a negative number.")
    } 
    if (init[2]<1){
      stop("Initial interval cannot be less than 1.")
    }
  } 
  zfit[1] <- init[1]
  xfit[1] <- init[2]
  
  if (length(w)==1){
    a.demand <- w[1]
    a.interval <- w[1]
  } else {
    a.demand <- w[1]
    a.interval <- w[2]
  }
  
  # Set coefficient
  if(type == "sba"){
    coeff <- 1-(a.interval/2)
  } else if(type == "sbj"){
    coeff <- (1-a.interval/(2-a.interval))
  } else {
    coeff <- 1
  }
  
  # Fit model
  for (i in 2:k){
    zfit[i] <- zfit[i-1] + a.demand * (z[i] - zfit[i-1])      # Demand
    xfit[i] <- xfit[i-1] + a.interval * (x[i] - xfit[i-1])    # Interval
  }
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
  
  # Prepare output - Assign weight to intervals if same with demand
  if (length(w)==1){
    w <- c(w,w)
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
  
  return(list(model=type,frc.in=frc.in,frc.out=frc.out,
              weights=w,initial=c(zfit[1],xfit[1]),components=comp))

}

#-------------------------------------------------
crost.opt <- function(data,type=c("croston","sba","sbj"),cost=c("mar","msr","mae","mse"),
                      w=NULL,nop=c(2,1),init,init.opt=c(TRUE,FALSE)){
# Optimisation function for Croston and variants
  
  type <- type[1]
  cost <- cost[1]
  nop <- nop[1]
  init.opt <- init.opt[1]
  
  # Croston decomposition
  nzd <- which(data != 0)               # Find location on non-zero demand
  k <- length(nzd)
  x <- c(nzd[1],nzd[2:k]-nzd[1:(k-1)])  # Intervals
  
  if (is.null(w) == TRUE && init.opt == FALSE){
    # Optimise only w
    p0 <- c(rep(0.05,nop))
    lbound <- c(rep(0,nop))
    ubound <- c(rep(1,nop))
    if (nop != 1){
      wopt <- optim(par=p0,crost.cost,method="Nelder-Mead",data=data,cost=cost,
                    type=type,w=w,nop=nop,w.opt=is.null(w),init=init,init.opt=init.opt,
                    lbound=lbound,ubound=ubound,control=list(maxit=2000))$par    
    } else {
      # Use Brent
      wopt <- optim(par=p0,crost.cost,method="Brent",data=data,cost=cost,
                    type=type,w=w,nop=nop,w.opt=is.null(w),init=init,init.opt=init.opt,
                    lbound=lbound,ubound=ubound,lower=lbound,upper=ubound,
                    control=list(maxit=2000))$par  
    }
    wopt <- c(wopt,init)
  } else if (is.null(w) == TRUE && init.opt == TRUE){
    # Optimise w and init
    p0 <- c(rep(0.05,nop),init[1],init[2])
    lbound <- c(rep(0,nop),0,1)
    ubound <- c(rep(1,nop),max(data),max(x))
    wopt <- optim(par=p0,crost.cost,method="Nelder-Mead",data=data,cost=cost,
                  type=type,w=w,nop=nop,w.opt=is.null(w),init=init,init.opt=init.opt,
                  lbound=lbound,ubound=ubound,control=list(maxit=2000))$par
  } else if (is.null(w) == FALSE && init.opt == TRUE){
    # Optimise only init
    nop <- length(w)
    p0 <- c(init[1],init[2])
    lbound <- c(0,1)
    ubound <- c(max(data),max(x))
    wopt <- optim(par=p0,crost.cost,method="Nelder-Mead",data=data,cost=cost,
                  type=type,w=w,nop=nop,w.opt=is.null(w),init=init,init.opt=init.opt,
                  lbound=lbound,ubound=ubound,control=list(maxit=2000))$par
    wopt <- c(w,wopt)
  }
  
  return(list(w=wopt[1:nop],init=wopt[(nop+1):(nop+2)]))
  
}

#-------------------------------------------------
crost.cost <- function(p0,data,cost,type,w,nop,w.opt,init,init.opt,lbound,ubound){
# Cost functions for Croston and variants  
  
  if (w.opt == TRUE && init.opt == TRUE){
    frc.in <- crost(data=data,w=p0[1:nop],h=0,init=p0[(nop+1):(nop+2)],
                    type=type,opt.on=TRUE)$frc.in
  } else if (w.opt == TRUE && init.opt == FALSE){
    frc.in <- crost(data=data,w=p0[1:nop],h=0,init=init,
                    type=type,opt.on=TRUE)$frc.in
  } else if (w.opt == FALSE && init.opt == TRUE){
    frc.in <- crost(data=data,w=w,h=0,init=p0,
                    type=type,opt.on=TRUE)$frc.in
  }
  
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
  
  # Constrains
  for (i in 1:(nop*w.opt+2*init.opt)){
    if (!(p0[i]>=lbound[i]) | !(p0[i]<=ubound[i])){
      E <- 9*10^99
    }
  }
  
  return(E)
  
}