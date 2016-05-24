###################################################
### chunk number 1: 
###################################################

algo.cusum <- function(disProgObj, control = list(range=range, k=1.04, h=2.26, m=NULL, trans="standard",alpha=NULL)){

  # Set the default values if not yet set
  if(is.null(control$k))
    control$k <- 1.04
  if(is.null(control$h))
    control$h <- 2.26
  if(is.null(control$trans))
    control$trans <- "standard"
    
  if(is.null(control$alpha))
    control$alpha <- 0.1
  alpha <- control$alpha

  observed <- disProgObj$observed
  timePoint <- control$range[1]

  # Estimate m (the expected number of cases), i.e. parameter lambda of a
  # poisson distribution based on time points 1:t-1
  if(is.null(control$m)) {
    m <- mean(observed[1:(timePoint-1)])
  } else if (is.numeric(control$m)) {
    m <- control$m
  } else if (control$m == "glm") {
    #Fit a glm to the first observations
    training <- 1:(timePoint-1)
    #Set the time index 
    t <- disProgObj$start[2] + training - 1
    #Set the observations
    x <- observed[training]
    #Set period
    p <- disProgObj$freq
    df <- data.frame(x=x,t=t)
    control$m.glm<- glm(x ~ 1 + cos(2*pi/p*t) + sin(2*pi/p*t) ,family=poisson(),data=df)

    #predict the values in range
    t.new <- disProgObj$start[2] + control$range - 1
    m <- predict(control$m.glm,newdata=data.frame(t=t.new),type="response")
  }


  #No transformation 
  #standObs <- observed[control$range]  
  x <- observed[control$range]
  standObs <- switch(control$trans,
    # compute standardized variables z3 (proposed by Rossi)
                     "rossi" =    (x - 3*m + 2*sqrt(x*m))/(2*sqrt(m)),
    # compute standardized variables z1 (based on asympotic normality)
                     "standard" = (x - m)/sqrt(m),
    # anscombe residuals
                     "anscombe" =  3/2*(x^(2/3)-m^(2/3))/m^(1/6),
                     
    # anscombe residuals as in pierce schafer based on 2nd order approx of E(X)
                     "anscombe2nd" =  (x^(2/3)-(m^(2/3)-m^(-1/3)/9))/(2/3*m^(1/6)),
                     
    # compute Pearson residuals for NegBin
                     "pearsonNegBin" = (x - m)/sqrt(m+alpha*m^2),
    # anscombe residuals for NegBin
                     "anscombeNegBin" =  anscombeNB(x,mu=m,alpha=alpha),
    # don't do anything
                     "none"     = x,
                     stop("invalid 'trans'formation")
                     )  

  # initialize the necessary vectors
  # start with cusum[timePoint -1] = 0, i.e. set cusum[1] = 0
  cusum <- matrix(0,nrow=(length(control$range)+1), ncol=1)
  alarm <- matrix(data = 0, nrow = (length(control$range)+1), ncol = 1)

  for (t in 1:length(control$range)){    
    # compute cumulated sums of standardized observations corrected with the
    # reference value k for all time points in range
    cusum[t+1]<- max(0, cusum[t]+(standObs[t]-control$k))
    
    # give alarm if the cusum is larger than the decision boundary h
    alarm[t+1] <- cusum[t+1] >= control$h
  }
  
  #Backtransform
  h <- control$h
  k <- control$k
  Ctm1 <-  cusum[1:length(control$range)]

  upperbound <- switch(control$trans,
    # standardized variables z3 (proposed by Rossi)
                     "rossi" =    2*h*m^(1/2)+2*k*m^(1/2)-2*Ctm1*m^(1/2)+5*m-2*(4*m^2+2*m^(3/2)*h+2*m^(3/2)*k-2*m^(3/2)*Ctm1)^(1/2),
    # standardized variables z1 (based on asympotic normality)
                     "standard" = ceiling(sqrt(m)*(h+k-Ctm1)+ m),
    # anscombe residuals
                     "anscombe" =  ifelse( ((2/3)*m^(1/6)*(h+k-Ctm1)+m^(2/3))<0, 
                                            0, 
                                           (2/3*m^(1/6)*(h+k-Ctm1)+m^(2/3))^(3/2) ),
    
    # anscombe residuals ?
                     "anscombe2nd" =  ifelse( ((2/3)*m^(1/6)*(h+k-Ctm1)+(m^(2/3)-m^(1/3)/9))<0, 
                                               0, 
                                              (2/3*m^(1/6)*(h+k-Ctm1)+(m^(2/3)-m^(1/3)/9))^(3/2) ),
    
    # Pearson residuals for NegBin  
                     "pearsonNegBin" = sqrt(m+alpha*m^2)*(h+k-Ctm1)+ m,
    # anscombe residuals for NegBin   ?
                     "anscombeNegBin" =  h-cusum[-1],
    # don't do anything
                     "none"     = h-cusum[-1]
                     )
  # ensure upper bound is positive and not NaN
  upperbound[is.na(upperbound)] <- 0
  upperbound[upperbound < 0] <- 0
    
  # discard cusum[1] and alarm[1]
  cusum <- cusum[-1]
  alarm <- alarm[-1]
  
  #Add name and data name to control object.
  control$name <- paste("cusum:", control$trans)
  control$data <- paste(deparse(substitute(disProgObj)))
  control$m    <- m

  # return alarm and upperbound vectors
  result <- list(alarm = alarm, upperbound = upperbound, disProgObj=disProgObj,control=control, cusum=cusum)

  class(result) = "survRes" # for surveillance system result
  return(result)
}


###################################################
### chunk number 2: 
###################################################
######################################################################
# Program to test the transformation of NegBin variables
# using the transformation similar to Anscombe residuals
######################################################################

#####################################################################
# function to evaluate hypgeom_2F1(1/3,2/3, 5/3, x)
# "exact" values for x = -(0:10) and linear interpolation for x = -(10:100)
####################################################################
hypgeom2F1special <- function(x) {
  #Return the z (the approximation grid), which is closest to x
  idx <- which.min(abs(surveillance.gvar.z-x))
  if(x >= -10)
    return(surveillance.gvar.hyp[idx])
  else{
    # find out interval that contains x
    if((x-surveillance.gvar.z[idx]) < 0){
      idxLow <- idx +1
      idxUp <- idx
    } else {
      idxLow <- idx
      idxUp <- idx -1
    }
    #linear interpolation: f(x)=f(x0)+(f(x1)-f(x0))/1*(x-x0)
    return(surveillance.gvar.hyp[idxLow]+(surveillance.gvar.hyp[idxUp]-surveillance.gvar.hyp[idxLow])*(x-surveillance.gvar.z[idxLow]))
  }
}

#####################################################################
# compute anscombe residuals for Y ~ NegBin(mu, alpha) using hypgeom2F1 function
# E(Y)= \mu, Var(Y) = \mu + \alpha*\mu^2
#################################################################
anscombeNB <- function(y,mu,alpha=0.1) {
  hypgeom.mu <- 3/2*mu^(2/3)*hypgeom2F1special(-alpha*mu)
  one <- function(y){
    up <- 3/2*y^(2/3) * hypgeom2F1special(-alpha*y) - hypgeom.mu
    down <- (mu+alpha*mu^2)^(1/6)
    return(up/down)
  }
  return(sapply(y,one))
}


###################################################
### chunk number 3: 
###################################################
######################################################################
# Given a specification of the average run length in the (a)cceptance
# and (r)ejected setting determine the k and h values in a standard
# normal setting.
#
# Description:
# Functions from the spc package are used in a simple univariate
# root finding problem.
#
# Params:
#  ARLa - average run length in acceptance setting (i.e. number before
#         false alarm
#  ARLw - average run length in rejection state (i.e. number before
#         an increase is detected (i.e. detection delay)
#  method - optim method to use, see ?optim
#
# Returns:
#  list( k - reference value, h - decision interval)
######################################################################
find.kh <- function(ARLa=500,ARLr=7,sided="one",method="BFGS",verbose=FALSE)
{
  if (!requireNamespace("spc"))
      stop("find.kh() requires package ", dQuote("spc"))
  #Small helper function which is to be minimized
  fun <- function(k) {
    if (k>0) {
      #Compute decision interval
      h <- spc::xcusum.crit(L0=ARLa,k=k,r=50,sided=sided)

      #Check if xcusum.crit managed to find a solution
      if (is.nan(h))
          stop("spc::xcusum.crit was not able to find a h corresponding to ",
               "ARLa=",ARLa," and k=",k)

      if (h > 0) {
        #Compute ARLr given the above computed h
        arlr <- spc::xcusum.arl(k,h,mu=2*k,r=50,sided=sided)
        #Deviation from the requested ARLr
        if (verbose) {
          cat("k=",k," score = ",(arlr-ARLr)^2,"\n")
        }
        return( (arlr-ARLr)^2 )
      } else {
        return(1e99)
      }
    } else {
      return( 1e99)
    }
  }
  k <- optim(1,fun,method=method)$par
  return(list(k=k,h=spc::xcusum.crit(L0=ARLa,k=k,r=50,sided=sided)))
}

