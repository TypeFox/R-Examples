#' Imbens-Kalyanaraman Optimal Bandwidth Calculation
#' 
#' \code{IKbandwidth} calculates the Imbens-Kalyanaraman optimal bandwidth
#' for local linear regression in Regression discontinuity designs.
#' 
#' @param X a numerical vector which is the running variable
#' @param Y a numerical vector which is the outcome variable
#' @param cutpoint the cutpoint
#' @param verbose logical flag indicating whether to print more information to the terminal. Default is \code{FALSE}.
#' @param kernel string indicating which kernel to use. Options are \code{"triangular"} (default and recommended), 
#' \code{"rectangular"}, \code{"epanechnikov"}, \code{"quartic"}, \code{"triweight"}, \code{"tricube"}, 
#' \code{"gaussian"}, and \code{"cosine"}.
#' @return The optimal bandwidth
#' @references Imbens, Guido and Karthik Kalyanaraman. (2009) "Optimal Bandwidth Choice for the regression discontinuity estimator," \emph{NBER Working Paper Series}. 14726. \url{http://www.nber.org/papers/w14726}
#' @importFrom stats complete.cases sd median lm coef
#' @export
#' @author Drew Dimmery <\email{drewd@@nyu.edu}>

IKbandwidth <-function (X,Y,cutpoint=NULL,verbose=FALSE,kernel="triangular") {
  #Implementation of Imbens-Kalyanaraman optimal bandwidth
  # for regression discontinuity
  sub<-complete.cases(X)&complete.cases(Y)
  X <- X[sub]
  Y <- Y[sub]
  Nx<-length(X)
  Ny<-length(Y)
  if(Nx!=Ny)
    stop("Running and outcome variable must be of equal length")
  if(is.null(cutpoint)) {
    cutpoint<-0
    if(verbose) cat("Using default cutpoint of zero.\n")
  } else {
    if(! (typeof(cutpoint) %in% c("integer","double")))
      stop("Cutpoint must be of a numeric type")
  }
  #Now we should be ready to start
  #Pilot bandwidth
  h1<-1.84*sd(X)*Nx^(-1/5)
  left<-X>=(cutpoint-h1) & X<=cutpoint
  right<-X>cutpoint & X<=(cutpoint+h1)
  Nl<-sum(left)
  Nr<-sum(right)
  Ybarl<-mean(Y[left])
  Ybarr<-mean(Y[right])
  fbarx<-(Nl+Nr)/(2*Nx*h1)
  varY<-(sum((Y[left]-Ybarl)^2)+sum((Y[right]-Ybarr)^2))/(Nl+Nr)
  medXl<-median(X[X<=cutpoint])
  medXr<-median(X[X>cutpoint])
  Nl<-sum(X<cutpoint)
  Nr<-sum(X>=cutpoint)
  cX<-X-cutpoint
  if(sum(X[left]>medXl)==0 | sum(X[right]<medXr)==0)
    stop("Insufficient data in vicinity of the cutpoint to calculate bandwidth.")
  #Model a cubic within the pilot bandwidth
  mod<-lm(Y~I(X>=cutpoint)+poly(cX,3,raw=T),subset=(X>=medXl&X<=medXr))
  m3<-6*coef(mod)[5]
  #New bandwidth estimate
  h2l<-3.56*(Nl^(-1/7))*(varY/(fbarx*max(m3^2,0.01)))^(1/7)
  h2r<-3.56*(Nr^(-1/7))*(varY/(fbarx*max(m3^2,0.01)))^(1/7)
  left<-(X>=(cutpoint-h2l)) & (X<cutpoint)
  right<-(X>=cutpoint) & (X<= (cutpoint+h2r))
  Nl<-sum(left)
  Nr<-sum(right)
  if(Nl==0 | Nr==0)
    stop("Insufficient data in vicinity of the cutpoint to calculate bandwidth.")
  #Estimate quadratics for curvature estimation
  mod<-lm(Y~poly(cX,2,raw=T),subset=right)
  m2r<-2*coef(mod)[3]
  mod<-lm(Y~poly(cX,2,raw=T),subset=left)
  m2l<-2*coef(mod)[3]
  rl<-720*varY/(Nl*(h2l^4))
  rr<-720*varY/(Nr*(h2r^4))
  #Which kernel are we using?
  # Method for finding these available in I--K p. 6
  if(kernel=="triangular") {
    ck<-3.43754
  } else if (kernel=="rectangular") {
    ck<-5.40384
  } else if(kernel=="epanechnikov") {
    ck<-3.1999
  } else if(kernel=="quartic" | kernel=="biweight") {
    ck<-3.65362
  } else if(kernel=="triweight") {
    ck<-4.06065
  } else if(kernel=="tricube") {
    ck<-3.68765
  } else if(kernel=="gaussian") {
    ck<-1.25864
  } else if(kernel=="cosine") {
    ck<-3.25869
  } else {
    stop("Unrecognized kernel.") 
  }
  #And there's our optimal bandwidth
  optbw<-ck*(2*varY/(fbarx*((m2r-m2l)^2+rr+rl)))^(1/5)*(Nx^(-1/5))
  left<-(X>=(cutpoint-optbw)) & (X<cutpoint)
  right<-(X>=cutpoint) & (X<= (cutpoint+optbw))
  if(sum(left)==0 | sum(right)==0)
    stop("Insufficient data in the calculated bandwidth.")
  names(optbw)<-NULL
  if(verbose) cat("Imbens-Kalyanamaran Optimal Bandwidth: ",sprintf("%.3f",optbw),"\n")
  return(optbw)
  }
