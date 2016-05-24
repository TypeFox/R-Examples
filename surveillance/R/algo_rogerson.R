###################################################
### chunk number 1: 
###################################################
###################################################################
# Average Run Lengths for CUSUMs using Markov Chain approach
#
# based on the program of Hawkins, D. M. (1992)
#   "Evaluation of Average Run Lengths of Cumulative Sum Charts
#    for an Arbitrary Data Distribution"
#   Communications in Statistics--Simulation. 21(4) 1001-1020.
#---------------------------------------------------------------
#
# for discrete distributions (i.e. Pois, Bin)
# and upward CUSUMS (increasing rate,probability)
#
# Parameters:
#  h - decision interval h
#  k - reference value k
#  distr - "poisson" or "binomial"
#  theta - distribution parameter for cdf distr, e.g. lambda for ppois, p for pbinomial
#  W - winsorizing value W (for robust CUSUM)
#      to get a nonrobust CUSUM set W > k+h
#  digits - k and h are rounded to digits decimal places
#  ... - further arguments for distribution
#        i.e. number of trials n for binomial (defaults to n=1)
#
# Returns:
#  ARL - one-sided ARL of the regular (no-head-start) CUSUM
###################################################################
arlCusum <- function(h=10, k=3, theta=2.4, distr=c("poisson","binomial"), W=NULL,digits=1,...){
  h <- round(h,digits)
  k <- round(k,digits)
  #cat("h",h,"k",k,"\n")

  distr <- match.arg(distr,c("poisson","binomial"))
  ##############
  # cdf of a binomial variate with fixed sample size
  pbinomial <- function(x,p,n=1){
    pbinom(x,size=n,prob=p)
  }
  ########
  distribution <- switch(distr,
                               "poisson" = ppois,
                               "binomial" = pbinomial
                        )

  #no winsorization
  if(is.null(W))
    W <- ceiling(h+k+1)

  # common denominator of h and k
  denrat <- commonDenom(h,k,digits=digits)
  #cat("h =",h,"k =",k,"denrat",denrat,"\n")

  # check parameters for feasibility
  if(h <=0)
    stop("Nonpositive decision interval\n")
  if(W <= k)
    stop("Winsorization value less than reference value\n")

  K <- as.integer(k*denrat+0.5)
  N <- as.integer(denrat)
  M <- as.integer(h*denrat -0.5)
  w <- as.integer(W*denrat+0.5)

  deviat <- abs(K-k*denrat)+abs(M-h*denrat+1)+abs(w-W*denrat)

  if(deviat > .01)
    stop("h, k or W not a multiple of 1/denrat\n")

  # determine probabilities
  x <- seq(((-M-1)+K)/N,(M+K)/N,by=(1/denrat))
  probs <- distribution(x, theta,...)

  # Winsorization (any observation exeeding some treshold value W is replaced by W
  # CUSUM is then: S_n = max{0, S_n-1 + min(X_n,W) - k}
  probs[x>=W] <- 1

  #construct transition matrix
  transition <- matrix(NA,M+1,M+1)
  transition[1,1] <- probs[(M+2)] #Pr(X <= k)
  transition[-1,1] <- probs[(M+2)-(1:M)]  #Pr(X <= -j+ k)  ,j=1,2,...,h-1
  transition[1,-1] <- probs[(M+2)+(1:M)]- probs[(M+2)+(1:M)-1] #Pr(X = j+ k)  , j=1,2,...,h-1

  idx <-rep((M+2):((M+2)+M-1),M)-rep(0:(M-1),each=M)
  transition[-1,-1] <- matrix(probs[idx]-probs[idx-1],nrow=M,ncol=M,byrow=TRUE)
  #print(transition)

  # I - transition matrix R
  IminusR <- diag(1,M+1) - transition

#Solve might work poorly in some cases
  res <- try(solve(IminusR)%*%rep(1,M+1),silent=TRUE)
#  res <- try(qr.solve(IminusR)%*%rep(1,M+1),silent=TRUE)
  if(inherits(res, "try-error")){
    warning("I-R singular\n")
    return(NA)
  }

  ARL <- res[1]
  #FIRARL - one-sided ARL of the FIR CUSUM with head start 0.5h
  FIRARL <- res[(M+1)/2+1]

  return(list(ARL=ARL,FIR.ARL=FIRARL))
}

#################################################################
# find smallest common denominator of x and y,
# i.e find an integer N so that x=X/N and y=Y/N (with X,Y integer)
#################################################################
commonDenom <- function(x,y,digits=1){
  x <- round(x,digits)
  y <- round(y,digits)

  which.max(  ( round((x+y)*1:(10^digits),digits)%%1 == 0 )     # (x+y)*N is integer
            & ( round(x*1:(10^digits),digits)%%1 == 0 )        # x*N is integer
            & ( round(y*1:(10^digits),digits)%%1 == 0 ) )      # y*N is integer
}



###################################################
### chunk number 2: 
###################################################
#################################################################
# find reference value k for a Poisson /Binomial CUSUM
# designed to detect a change from theta0 to theta1
#
# digits - k is rounded to digits decimal places if roundK=TRUE
# ... - extra arguments for distribution,
#       i.e number of trials n for binomial, set to 1 if not specified
##################################################################
findK <- function(theta0,theta1,distr=c("poisson","binomial"),roundK=FALSE,digits=1,...){
  n <- list(...)$n
  if(is.null(n))
    n <- 1
  distr <- match.arg(distr,c("poisson","binomial"))
  k <- switch(distr,
         "poisson" = (theta1 - theta0)/(log(theta1)-log(theta0)),
         "binomial" = -n*(log(1-theta1)-log(1-theta0))/(log(theta1*(1-theta0))-log(theta0*(1-theta1)))
         )

  # for discrete data the
  # Cusum values are of form integer - integer multiple of k
  # so there is a limited set of possible values of h (and ARL)

  if(roundK){
  # add/substract 0.05 to k so that k isn't an integer or a multiple of 0.5
  # when rounded (then the Markov Chain has more states)
    if(round(k,digits=digits)%%1 == 0.5 | round(k,digits=digits)%%1 == 0){
      round.k <- ((k-floor(k))*10^digits)%%1
      #print(roundK)
        if(round.k < .5 )
          k <- k+0.5*10^(-digits)
        else
          k <- k-0.5*10^(-digits)
    }
    k <- round(k,digits=digits)
  }
  return(k)
}


###################################################
### chunk number 3: 
###################################################
##################################################################
# function to find the decision limit h so that the
# average run length for a Poisson/Binomial CUSUM with in-control
# parameter theta0 is (approx.) ARL0
#
# Params:
#  ARL0 - desired in-control ARL
#  theta0 - in-control parameter
#  s - change to detect (in stdev)
#  rel.tol - (relative) tolerance (if attainable)
#  roundK - should k be rounded up to digits decimal place (avoiding integers, multiples of 0.5)
#  digits - h is rounded to digits decimal places
#  distr - "poisson" or "binomial"
#  ... - further arguments for distribution (i.e number of trials n for "binomial")
#
# Returns:
#  vector c(theta0, h, k, ARL, rel.tol)
#################################################################
findH <- function(ARL0,theta0,s=1, rel.tol=0.03,roundK=TRUE,distr=c("poisson","binomial"),digits=1,FIR=FALSE,...){
  distr <- match.arg(distr,c("poisson","binomial"))
  #FIR-ARL or zero-start ARL?
  fir.arl <- ifelse(FIR,2,1)

  theta1 <- getTheta1(theta0,s=s,distr=distr)

  k <- findK(theta0,theta1,roundK=roundK,distr=distr,digits=digits,...)

  # compute ARL for two (arbitrary) points (h1, h2)
  h1 <- min(12,4*k)
  arl1 <- arlCusum(h=h1,k=k,theta=theta0,distr=distr,digits=digits,...)[[fir.arl]]
  nEval <- 1

  #ensure h1 and arl1 are not too small (-> log. interpolation is better)
  while(arl1 < 100){
    h1 <- 2*h1
    arl1 <- arlCusum(h=h1,k=k,theta=theta0,distr=distr,digits=digits,...)[[fir.arl]]
    nEval <- nEval + 1
  }
  h2 <- h1*2^(sign(ARL0-arl1))
  arl2 <- arlCusum(h=h2,k=k,theta=theta0,distr=distr,digits=digits,...)[[fir.arl]]
  nEval <- nEval + 1

  # determine h (that leads to an ARL of ARL0) using logarithmic interpolation
  h.hat <- round(logInterpolation(ARL0,h1,h2,arl1,arl2),digits)

  # what's the actual ARL for h
  arl <- arlCusum(h=h.hat,k=k,theta=theta0,distr=distr,digits=digits,...)[[fir.arl]]
  nEval <- nEval + 1

  relTol <- abs((arl-ARL0)/ARL0)
  #cat("theta0:", theta0,"k:", k,"h:", h.hat,"ARL:",arl,"relTol:", relTol,"\n")

  i<-0
  signs <- sign(ARL0-arl)
  convergence <- relTol < rel.tol
  if(convergence){
#    print(nEval)
    return(c("theta0"=theta0,"h"=h.hat,"k"=k,"ARL"=arl,"rel.tol"=relTol))
  }

  # find hLow so that the target ARL0 is in interval c(ARL(hLow), ARL(h.hat))
  denrat <- 1/commonDenom(1,k,digits=digits)
  steps <- denrat #max(0.1,denrat)
#  cat("denrat",denrat,"steps",steps,"\n")

  hLow <- round(h.hat+signs*steps,digits)
  arlLow <- arlCusum(h=hLow,k=k,theta=theta0,distr=distr,digits=digits,...)[[fir.arl]]
  nEval <- nEval + 1
  relTol.Low <- abs((arlLow-ARL0)/ARL0)
  if(relTol.Low < rel.tol){
#    print(nEval)
    return(c("theta0"=theta0,"h"=hLow,"k"=k,"ARL"=arlLow,"rel.tol"=relTol.Low))
  }

  while(sign(ARL0-arl)*sign(ARL0-arlLow)>0){
#    cat("steps:",nEval,"\n")
    h.hat <- hLow
    arl <-arlLow
    relTol <- relTol.Low
    signs <- sign(ARL0-arl)

    hLow <- round(h.hat+signs*steps,digits)
    arlLow <- arlCusum(h=hLow,k=k,theta=theta0,distr=distr,digits=digits,...)[[fir.arl]]
    nEval <- nEval + 1
    relTol.Low <- abs((arlLow-ARL0)/ARL0)
    if(relTol.Low < rel.tol){
#      print(nEval)
      return(c("theta0"=theta0,"h"=hLow,"k"=k,"ARL"=arlLow,"rel.tol"=relTol.Low))
    }
#    cat("hLow:", hLow,"ARL:",arlLow,"relTol:",relTol.Low,"\n")
  }
#  cat("hLow:", hLow,"ARL:",arlLow,"relTol:",relTol.Low,"\n")

  # return the ARL which is at least the target ARL0
  if(sign(ARL0-arlLow)<0){
    h.hat <- hLow
    arl <- arlLow
    relTol <- relTol.Low
  }
  #print(nEval)
  return(c("theta0"=theta0,"h"=h.hat,"k"=k,"ARL"=arl,"rel.tol"=relTol))
}

##################################################################
# find h for various values theta0
#
# Params:
#  theta0 - vector of in control parameter
#  ARL0 - desired in-control ARL
#
# Returns:
# matrix with columns c(theta0, h, k, ARL, rel.Tol)
##################################################################
hValues <- function(theta0,ARL0,rel.tol=0.02,s=1,roundK=TRUE,digits=1,distr=c("poisson","binomial"),FIR=FALSE,...){
  distr <- match.arg(distr,c("poisson","binomial"))
  n <- list(...)$n
  hVals <- t(sapply(theta0,findH,ARL0=ARL0,rel.tol=rel.tol,s=s,roundK=roundK,digits=digits,distr=distr,FIR=FIR,...))
  res <- list(hValues=hVals,ARL0=ARL0,s=s,rel.tol=rel.tol,distribution=distr,firARL=FIR)
  res$n <- n
  return(res)
}



##################################################################
# get the decision limit h for CUSUM with
# in-control parameter theta using a "table" of h values
#
#  theta - in-control parameter
#  hValues - matrix with columns c(theta, h)
##################################################################
getH <- function(theta,hValues){
  one<- function(theta){
    theta.diff <- abs(hValues[,1]-theta)
    idx <- which.min(theta.diff)
    h <- hValues[idx,2]

    if(theta.diff[idx] > 0.05)
      warning("table doesn't contain h value for theta = ",theta,"\n")

    return(h)
  }
  sapply(theta,one)
}

#####################################################################
# get decision interval h and reference value k
#####################################################################
getHK <- function(theta,hValues){
  one<- function(theta){
    theta.diff <- abs(hValues[,1]-theta)
    idx <- which.min(theta.diff)
    hk <- hValues[idx,2:3]

    if(theta.diff[idx] > 0.05)
      warning("table doesn't contain h value for theta = ",theta,"\n")

    return(hk)
  }
  t(sapply(theta,one))
}


#################################################################
# get out-of-control parameter theta1
#
# X ~ Po(lambda0): theta1 = lambda0 + s*sqrt(lambda0)
#  theta1 corresponds to a s*stdev increase  in mean
#
# X ~Bin(n,pi)
#  H0: Odds of failure =pi/(1-pi)  vs  H1: Odds = s*pi/(1-pi)
# prob of failure under H1 is then pi1 = s*pi0/(1+(s-1)*pi0)
#################################################################
getTheta1 <- function(theta0,s=1,distr=c("poisson","binomial")){
  distr <- match.arg(distr,c("poisson","binomial"))
  theta1 <- switch(distr,
                   "poisson" = theta0 + s*sqrt(theta0),
                   "binomial" = s*theta0/(1-theta0+s*theta0)
                   )
  return(theta1)
}




#################################################################
# logarithmic interpolation, i.e. linear interpolation of ln(f(x))
# in points (x0,ln(f0)), (x1,ln(f1))
#
#   (ln(f)-ln(f0))/(ln(f1)-ln(f0)) = (x-x0)/(x1-x0)
#
# returns: x
#
# to find decision limit h for given ARL0 set x = h, f(h) = ARL0(h,k)
# and solve equation for x
#################################################################
logInterpolation <- function(f,x0,x1,f0,f1){
  x0 + ((x1-x0)*(log(f)-log(f0)))/(log(f1)-log(f0))
}



###################################################
### chunk number 4: 
###################################################
#  control - list with
#     range - vector of indices in the observed matrix to monitor
#     theta0t - matrix with in-control parameter, needs to be specified
#     ARL0 - desired average run length for each one of the univariate CUSUMs
#     s -  change to detect
#     hValues - matrix with decision intervals for theta0_t
#     reset - if TRUE, the CUSUM is reset to zero after an alarm
#     nt - time varying sample sizes (for Binomial),
#          matrix of same dimension as theta0t
algo.rogerson <- function(disProgObj,
             control=list(range=range, theta0t=NULL, ARL0=NULL, s=NULL, hValues=NULL,
             distribution=c("poisson","binomial"), nt=NULL, FIR=FALSE,limit=NULL, digits=1)){

  if (is.null(control$s)) { stop("Error: the s value is not specified") }
  if (is.null(control$hValues)) { stop("Error: the hValues are not specified") }
#  if (is.null(control$ARL0)) { stop("Error: no ARL0 value specified") }

  #Default value is poisson
  control$distribution <- match.arg(control$distribution,c("poisson","binomial"))
  if(is.null(control$FIR)){
    control$FIR <- FALSE
  }
  if(is.null(control$limit))
    control$limit <- -1
  if(is.null(control$digits))
    control$digits <- 1

  x <- as.matrix(disProgObj$observed[control$range,])

  if (is.null(control$theta0t)) {
    stop("Error: no theta0t vector specified")
  } else {
    theta0t <- as.matrix(control$theta0t)
  }
  #theta0 <- colMeans(theta0t)

  #size = length of process
  size <- nrow(x)
  nAreas <- ncol(theta0t)
  theta0 <- rep(mean(theta0t),nAreas)


  #check dimensions of x, theta0t
  if(size !=nrow(theta0t) | (ncol(x)%%nAreas)!=0)
    stop("wrong dimensions\n")

  reps <- ncol(x)/nAreas

  #time-varying size n for Binomial
  nt<-control$nt
  if(control$distribution=="binomial"){
    if(is.null(nt))
      nt <- matrix(rep(control$n,size),ncol=1)
    else
      nt<-as.matrix(nt)
  }

  theta1 <- getTheta1(theta0,s=control$s,distr=control$distribution)
  theta1t <- getTheta1(theta0t,s=control$s,distr=control$distribution)

  hk <- getHK(theta0,hValues=control$hValues)
  k <- hk[,"k"]
  h <- hk[,"h"]
  #cat("k =",k,"h =",h,"\n")

  if(control$FIR){
    control$limit <- 0.5
    fir <- h/2
  } else {
    fir <- 0
  }
  #cat("fir",fir,"\n")

  # initialize the necessary vectors
  # start with cusum[1] = 0
  cusum <- matrix(0,nrow=(size+1), ncol=nAreas*reps)
  cusum[1,] <- fir
  alarm <- matrix(data = 0, nrow = (size+1), ncol = nAreas*reps)
  upperbound <- matrix(0,nrow=(size+1),ncol=reps)

  #CUSUM as in Rogerson (2004)
  for(t in 1:size){
    #choose k_t based upon theta_0t and theta_1t
    hkt <- getHK(theta0t[t,],hValues=control$hValues)
    #kt <- hkt[,"k"]
    kt <- findK(theta0t[t,],theta1t[t,],distr=control$distribution,roundK=TRUE,
                digits=control$digits, n=nt[t,])   #

    #for given k_t (theta0t) and ARL_0 choose h_t
    #ht <- getH(lambda0t[t],control$hValues)
    ht <- hkt[,"h"]
    ct <- h/ht

    # compute cumulative sums of observations x corrected with the
    # reference value kt, scaled by factor ct
#    cusum[t+1,]<- pmax(0, cusum[t,] + ct*(x[t,]-kt))

    # reset CUSUM to zero if an alarm is given at time t
    if((control$limit >= 0) & any(alarm[t,]==1)){
      cusum.t <- cusum[t,]
      cusum.t[alarm[t,]==1] <- pmin(cusum[t,], control$limit*h)[alarm[t,]==1]
      cusum[t+1,]<- pmax(0, cusum.t + ct*(x[t,]-kt))
    } else {
      cusum[t+1,]<- pmax(0, cusum[t,] + ct*(x[t,]-kt))
    }
    # give alarm if the cusum is larger than h
    alarm[t+1,] <- cusum[t+1,] >= h
    # in case speed is premium then one might want to comment this line
    if((control$limit >= 0) & any(alarm[t,]==1)) {
      upperbound[t+1,] <- ceiling( (h-cusum.t)/ct + kt) 
    } else {
      upperbound[t+1,] <- ceiling( (h-cusum[t,])/ct + kt) 
    }

    #Ensure upperbound is positive (this should always be the case)
    if (upperbound[t+1,] < 0) { upperbound[t+1,] <- 0}    
  }


  # discard cusum[1] and alarm[1]
  cusum <- as.matrix(cusum[-1,])
  alarm <- as.matrix(alarm[-1,])
  upperbound <- as.matrix(upperbound[-1,])

  #Add name and data name to control object.
  control$name <- paste("CUSUM Rogerson:",control$distribution)
  control$data <- paste(deparse(substitute(disProgObj)))

  # return alarm and upperbound vectors
  result <- list(alarm = alarm, upperbound = upperbound, disProgObj=disProgObj,control=c(control,list(h=h)))
  class(result) = "survRes" # for surveillance system result
  return(result)
}






