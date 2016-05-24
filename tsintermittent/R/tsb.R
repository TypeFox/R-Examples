tsb <- function(data,h=10,w=NULL,init=c("mean","naive"),
                cost=c("mar","msr","mae","mse"),
                init.opt=c(TRUE,FALSE),outplot=c(FALSE,TRUE),
                opt.on=c(FALSE,TRUE),na.rm=c(FALSE,TRUE)){
# TSB method
#
# Inputs:
#   data        Intermittent demand time series.
#   h           Forecast horizon.
#   w           Smoothing parameters. If w == NULL then parameters are optimised.
#               Otherwise first parameter is for demand and second for demand probability. 
#   init        Initial values for demand and intervals. This can be:
#                 c(z,x)  - Vector of two scalars, where first is initial demand and 
#                           second is initial interval;
#                 "naive" - Initial demand is first non-zero demand and initial demand
#                           probability is again the first one;
#                 "mean"  - Same as "naive", but initial demand probability is the mean 
#                           of all in sample probabilities.
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
#   weights     Smoothing parameters for demand and demand probability.
#   initial     Initialisation values for demand and demand probability smoothing.
#
# Example:
#   tsb(ts.data1,outplot=TRUE)
#
# Notes:
# Optimisation of the method described in:
# N. Kourentzes, 2014, On intermittent demand model optimisation and selection, 
# International Journal of Production Economics, 156: 180-190. 
# http://dx.doi.org/10.1016/j.ijpe.2014.06.007
# http://kourentzes.com/forecasting/2014/06/11/on-intermittent-demand-model-optimisation-and-selection/
#  
# Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>
  
  # Defaults
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
  
  # TSB decomposition
  p <- as.numeric(data!=0)    # Demand probability
  z <- data[data!=0]          # Non-zero demand
  
  # Initialise
  if (!(is.numeric(init) && length(init)==2)){
    if (init=="mean"){
      init <- c(z[1],mean(p))
    } else {
      init <- c(z[1],p[1])
    }
  }
  
  # Optimise parameters if requested
  if (opt.on == FALSE){
    if (is.null(w) || init.opt == TRUE){
      wopt <- tsb.opt(data,cost,w,init,init.opt)
      w <- wopt$w
      init <- wopt$init
    } else {
      if (length(w)!=2){
        stop(paste("w must be a vector of 2 elements: the smoothing parameter",
                   " for the non-zero demand and the parameter for the ",
                   "probability of demand.",sep=""))
      }
    }
  }
  
  # Pre-allocate memory
  zfit <- vector("numeric",n)
  pfit <- vector("numeric",n)
  
  # Assign initial values and parameters
  if (opt.on == FALSE){
    if (init[1]<0){
      stop("Initial demand cannot be a negative number.")
    } 
    if (init[2]<0){
      stop("Initial demand probability cannot be a negative number.")
    }
  }
  zfit[1] <- init[1]
  pfit[1] <- init[2]
  
  # Fit model
  for (i in 2:n){
    pfit[i] <- pfit[i-1] + w[2]*(p[i]-pfit[i-1])        # Demand probability
    if (p[i]==0){
      zfit[i] <- zfit[i-1]
    } else {
      zfit[i] <- zfit[i-1] + w[1]*(data[i]-zfit[i-1])   # Demand
    }
  }
  yfit <- pfit*zfit
  
  frc.in <- c(NA,yfit[1:(n-1)])
  if (h>0){
    frc.out <- rep(yfit[n],h)
  } else {
    frc.out = NULL
  }
  
  # Plot
  if (outplot==TRUE){
    plot(1:n,data,type="l",xlim=c(1,(n+h)),xlab="Period",ylab="",
         xaxs="i",yaxs="i",ylim=c(0,max(data)*1.1))
    lines(which(data>0),data[data>0],type="p",pch=20)
    lines(1:n,frc.in,col="red")
    lines((n+1):(n+h),frc.out,col="red",lwd=2)
  }
  
  return(list(model="tsb",frc.in=frc.in,frc.out=frc.out,
              weights=w,initial=c(zfit[1],pfit[1])))
  
}

#-------------------------------------------------
tsb.opt <- function(data,cost=c("mar","msr","mae","mse"),w=NULL,
                    init,init.opt=c(TRUE,FALSE)){
# Optimisation function for TSB
  
  cost <- cost[1]
  init.opt <- init.opt[1]
  
  if (is.null(w) == TRUE && init.opt == TRUE){
    # Optimise w and init
    p0 <- c(rep(0.05,2),init[1],init[2])
    lbound <- c(0,0,0,0)
    ubound <- c(1,1,max(data),1)
    wopt <- optim(par=p0,tsb.cost,method="Nelder-Mead",data=data,cost=cost,
                  w=w,w.opt=is.null(w),init=init,init.opt=init.opt,
                  lbound=lbound,ubound=ubound,control=list(maxit = 2000))$par
  } else if (is.null(w) == TRUE && init.opt == FALSE){
    # Optimise only w
    p0 <- c(rep(0.05,2))
    lbound <- c(0,0)
    ubound <- c(1,1)
    wopt <- optim(par=p0,tsb.cost,method="Nelder-Mead",data=data,cost=cost,
                  w=w,w.opt=is.null(w),init=init,init.opt=init.opt,
                  lbound=lbound,ubound=ubound,control=list(maxit = 2000))$par   
    wopt <- c(wopt,init)
  } else if (is.null(w) == FALSE && init.opt == TRUE){
    # Optimise only init
    p0 <- c(init[1],init[2])
    lbound <- c(0,0)
    ubound <- c(max(data),1)
    wopt <- optim(par=p0,tsb.cost,method="Nelder-Mead",data=data,cost=cost,
                  w=w,w.opt=is.null(w),init=init,init.opt=init.opt,
                  lbound=lbound,ubound=ubound,control=list(maxit = 2000))$par
    wopt <- c(w,wopt)
  }
  
  return(list(w=wopt[1:2],init=wopt[3:4]))
  
}

#-------------------------------------------------
tsb.cost <- function(p0,data,cost,w,w.opt,init,init.opt,lbound,ubound){
  # Cost functions for TSB
  
  if (w.opt == TRUE && init.opt == TRUE){
    frc.in <- tsb(data=data,w=p0[1:2],h=0,init=p0[3:4],opt.on=TRUE)$frc.in
  } else if (w.opt == TRUE && init.opt == FALSE){
    frc.in <- tsb(data=data,w=p0[1:2],h=0,init=init,opt.on=TRUE)$frc.in
  } else if (w.opt == FALSE && init.opt == TRUE){
    frc.in <- tsb(data=data,w=w,h=0,init=p0[1:2],opt.on=TRUE)$frc.in
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
  
  # Bounds
  for (i in 1:(2*w.opt+2*init.opt)){
    if (!p0[i]>=lbound[i] | !p0[i]<=ubound[i]){
      E <- 9*10^99
    }
  }
  
  # Parameter of demand probability must be smaller than parameter of demand
  if (w.opt == TRUE){
    if (p0[1] < p0[2]){
      E <- 9*10^99
    }
  }
  
  return(E)
  
}