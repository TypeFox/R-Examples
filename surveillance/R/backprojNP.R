######################################################################
# Implementation of the backprojection method as described in
# Becker et al. (1991), Stats in Med, 10, 1527-1542. The method
# was originally developed for the back-projection of AIDS incidence
# but it is equally useful for analysing the epidemic curve in outbreak
# situations of a disease with long incubation time, e.g. in order
# to illustrate the effect of intervention measures.
#
# See backprojNP.Rd for the remaining details.
######################################################################

######################################################################
# Helper function: Replace NaN or is.infinite values with zero.
# Good against division by zero problems.
#
# Parameters:
#  x - a vector of type double
######################################################################
naninf2zero <- function(x) {x[is.nan(x) | is.infinite(x)] <- 0 ; return(x)}


## Rcpp inline function to significantly speed up the computation of equation
## 3a in the Becker et al. (1991) paper. Created with the help of Daniel
## Sabanes Bove, University of Zurich.

## eq3a <-
##     cxxfunction(signature(rlambdaOld="numeric",
##                           ry="numeric",
##                           rincuPmf="numeric"),
##                 '
## {
##     // get arguments
## NumericVector lambdaOld(rlambdaOld);
## int T = lambdaOld.length();
## NumericVector y(ry);
## NumericVector incuPmf(rincuPmf);
## NumericVector dincu(T);
## NumericVector pincu(T);
## pincu[0] = dincu[0];    
## for (int i=1; i<incuPmf.length(); i++) {
##   dincu[i] = incuPmf[i];
##   pincu[i] = pincu[i-1] + dincu[i];
## }
## for (int i=incuPmf.length(); i<T; i++) {
##   dincu[i] = 0.0;
##   pincu[i] = 1.0;
## }  

      
## // result vector
## NumericVector phiNew(T);

## // many loops
## for(int t = 0; t < T; ++t)
## {
##     double sum3a = 0.0;

##     for(int d = 0; d <= T - (t + 1); ++d)
##     {
## 	double tmp = 0.0;

## 	for(int i = 0; i < t + d; ++i)
## 	{
## 	    tmp += lambdaOld[i] * dincu[t + d - i];
## 	}
       
## 	tmp = dincu[d] / tmp;

## 	if(R_IsNaN(tmp) || (! R_finite(tmp)))
## 	{
## 	    tmp = 0.0;
## 	}

## 	sum3a += y[t + d] * tmp;	    
##     }
    

##     //printf("Querying index %d\\n", T-(t+1));
    
##     phiNew[t] = lambdaOld[t] / pincu[T - (t + 1)] * sum3a;
    
##     if(R_IsNaN(phiNew[t]) || (! R_finite(phiNew[t])))
##     {
## 	phiNew[t] = 0.0;
##     }
## }

##     //Show values
##     //for (int i=0; i<T; ++i) {
##     //  printf( "%f\\n", phiNew[i]);
##     //}

## return wrap(phiNew);} ',
## plugin="Rcpp")


######################################################################
# Single step of the EMS algorithm by Becker et al (1991). This function
# is called by backprojNP.
#
# Parameters:
#  lambda.old - vector of length T containing the current rates
#  Y          - vector of length T containing the observed values
#  dincu      - probability mass function of the incubation time. I.e.
#               a function to be evaluated at integer numbers
#  pincu      - cumulative mass function of the incubation time, i.e. an
#               object of type function. Needs to in sync with dincu. 
#  k          - smoothing parameter of the EMS algo,
#               needs to be an even number
#
# Returns:
# 
######################################################################

em.step.becker <- function(lambda.old, Y, dincu, pincu, k, incu.pmf, eq3a.method=c("R","C")) {
  #k needs to be divisible by two
  if (k %% 2 != 0) stop("k needs to be even.")
  #which method to use 
  eq3a.method <- match.arg(eq3a.method,c("R","C"))

  #Initialize
  T <- length(Y)
  
  #Define new parameters
  phi.new <- lambda.new <- 0*lambda.old

  if (eq3a.method=="R") {
    #EM step. Problem that some of the sums can be zero if the incubation
    #distribution has zeroes at d=0,1,2
     for (t in 1:T) {
    #Calculate sum as in equation (3a) of Becker (1991)
      sum3a <- 0
      for (d in 0:(T-t)) {
        sum3a <- sum3a + Y[t+d] * naninf2zero(dincu(d) / sum(sapply(1:(t+d),function(i) lambda.old[i]*dincu(t+d-i))))
      }
      phi.new[t] <- naninf2zero(lambda.old[t]/pincu(T-t)) * sum3a
    }
  } else {
    phi.new <- .Call("eq3a",lambda.old=as.numeric(lambda.old),Y=as.numeric(Y),incu.pmf=as.numeric(incu.pmf),PACKAGE="surveillance")
    #phi.new <- eq3a(lambda.old, Y, incu.pmf)
  }

  
  #Smoothing step
  if (k>0) {
    w <- choose(k,0:k)/2^k
    for (t in 1:T) {
      i.sub <- t+(0:k)-k/2
      goodIdx <- i.sub %in% 1:T
      w.sub <- w[goodIdx]/sum(w[goodIdx])
      lambda.new[t] <- sum(w.sub * phi.new[i.sub[goodIdx]])
    }
  } else { #no smoothing
    lambda.new <- phi.new
  }
  
  #Done.
  return(lambda=lambda.new)
}

######################################################################
# STS compatible function to call the non-parametric back-projection
# method of Becker et al (1991) for time aggregated data.
#
# Parameters:
#  sts - sts object with the observed incidence as "observed" slot
#  incu.pmf - incubation time pmf as a vector with index 0,..,d_max. Please
#                 note that the support includes zero!
#  k - smoothing parameter for the EMS algorithm
#  eps - relative convergence criteration
#  iter.max - max number of iterations
#  verbose - boolean, if TRUE provide extra output when running the method
#  lambda0 - start value for lambda, default: uniform
#  hookFun - hook function to call after each EMS step, a function
#            of type hookFun=function(stsj,...)
#
# Returns:
#  sts object with upperbound set to the backprojected lambda.
######################################################################


backprojNP.fit <- function(sts, incu.pmf,k=2,eps=1e-5,iter.max=250,verbose=FALSE,lambda0=NULL,eq3a.method=c("R","C"),hookFun=function(stsbp) {}, ...) {

  #Determine method
  eq3a.method <- match.arg(eq3a.method, c("R","C"))

  #Define object to return
  lambda.hat <- matrix(NA,ncol=ncol(sts),nrow=nrow(sts))

  #Loop over all series
  for (j in 1:ncol(sts)) {
    #Inform (if requested) what series we are looking at
    if ((ncol(sts)>1) & verbose) {
      cat("Backprojecting series no. ",j,"\n")
    }

    #Extract incidence time series
    Y <- observed(sts)[,j]

    #If default behaviour for lambda0 is desired
    if (is.null(lambda0)) {
      lambda0j <- rep(sum(Y)/length(Y),length(Y))
    } else {
      lambda0j <- lambda0[,j]
    }

    #Create incubation time distribution vectors for the j'th series
    inc.pmf <- as.numeric(incu.pmf[,j])
    inc.cdf <- cumsum(inc.pmf)
  
    #Create wrapper functions for the PMF and CDF based on the vector.
    #These function will be used in the R version of eq3a.
    #ToDo: The function uses the global variable inc.pmf which
    #definitely is dirty coding. How to define this function
    #in an environment where inc.pmf is present?
    dincu <- function(x) {
      notInSupport <- x<0 | x>=length(inc.pmf)
      #Give index -1 to invalid queries
      x[notInSupport] <- -1
      return(c(0,inc.pmf)[x+2])
    }
    #Cumulative distribution function. Uses global var "inc.cdf"
    pincu <- function(x) {
      x[x<0] <- -1
      x[x>=length(inc.cdf)] <- length(inc.cdf)-1
      return(c(0,inc.cdf)[x+2])
    }
   
    #Iteration counter and convergence indicator
    i <- 0
    stop <- FALSE
    lambda <- lambda0j
  
    #Loop until stop 
    while (!stop) {
      #Add to counter
      i <- i+1
      lambda.i <- lambda
      #Perform one step
      lambda <- em.step.becker(lambda.old=lambda.i,Y=Y,dincu=dincu,pincu=pincu,k=k, incu.pmf=inc.pmf, eq3a.method=eq3a.method)
      
      #check stop
      #In original paper the expression to do so appears funny since
      #- and + deviations cancel. More realistic:
      #criterion <- abs(sum(res$lambda) - sum(lambda.i))/sum(lambda.i)
      criterion <- sqrt(sum((lambda- lambda.i)^2))/sqrt(sum(lambda.i^2))

      if (verbose) {
        cat("Convergence criterion @ iteration i=",i,": ", criterion,"\n")
      }
      #Check whether to stop
      stop <- criterion < eps | (i>iter.max)
      
      #Call hook function
      stsj <- sts[,j]
      upperbound(stsj) <- matrix(lambda,ncol=1)
      hookFun(stsj, ...)
    }
    #Done
    lambda.hat[,j] <- lambda
  }

  #Create new object with return put in the lambda slot
  bp.sts <- as(sts,"stsBP")
  bp.sts@upperbound <- lambda.hat
  bp.sts@control <- list(k=k,eps=eps,iter=i)
  return(bp.sts)
}


######################################################################
# EMS back-projection method including bootstrap based confidence
# intervals. The theory is indirectly given in Becker and Marschner (1993),
# Biometrika, 80(1):165-178 and more specifically in Yip et al, 2011,
# Communications in Statistics -- Simulation and Computation,
# 37(2):425-433.
#
# Parameters:
#
#  sts - sts object with the observed incidence as "observed" slot
#  incu.pmf - incubation time pmf as a vector with index 0,..,d_max. Please
#                 note that the support includes zero!
#  k - smoothing parameter for the EMS algorithm
#  eps - relative convergence criteration. If a vector of length two
#        then the first argument is used for the k=0 initial fit and
#        the second element for all EMS fits
#
#  iter.max - max number of iterations. Can be a vector of length two.
#             Similar use as in eps.
#  verbose - boolean, if TRUE provide extra output when running the method
#  lambda0 - start value for lambda, default: uniform
#  hookFun - hook function to call after each EMS step, a function
#            of type hookFun=function(Y,lambda,...)
#  B - number of bootstrap replicates. If B=-1 then no bootstrap CIs
#      are calculated.
#
# Returns:
#  sts object with upperbound set to the backprojected lambda.
######################################################################

backprojNP <- function(sts, incu.pmf,control=list(k=2,eps=rep(0.005,2),iter.max=rep(250,2),Tmark=nrow(sts),B=-1,alpha=0.05,verbose=FALSE,lambda0=NULL,eq3a.method=c("R","C"),hookFun=function(stsbp) {}),...) {

  #Check if backprojection is to be done multivariate time series case.
  if (ncol(sts)>1) {
    warning("Multivariate time series: Backprojection uses same eps for the individual time series.")
  }
  #Check if incu.pmf vector fits the dimension of the sts object. If not
  #either replicate it or throw an error.
  if (is.matrix(incu.pmf)) {
    if (!ncol(incu.pmf) == ncol(sts)) {
      stop("Dimensions of sts object and incu.pmf don't match.")
    }
  } else {
    if (ncol(sts)>1) {
      warning("Backprojection uses same incubation time distribution for the individual time series.")
    }
    incu.pmf <- matrix(incu.pmf,ncol=ncol(sts),dimnames=list(NULL,colnames(sts)))
  }

  #Fill control object as appropriate and in sync with the default value
  if (is.null(control[["k",exact=TRUE]]))  { control$k <- 2 }
  if (is.null(control[["eps",exact=TRUE]])) { control$eps <- rep(0.005,2) }
  if (is.null(control[["iter.max",exact=TRUE]])) { control$iter.max <- rep(250,2) }
  if (is.null(control[["Tmark",exact=TRUE]])) { control$Tmark <- nrow(sts) }
  if (is.null(control[["B",exact=TRUE]])) { control$B <- -1 }
  if (is.null(control[["alpha",exact=TRUE]])) { control$alpha <- 0.05 }
  if (is.null(control[["verbose",exact=TRUE]])) { control$verbose <- FALSE }
  if (is.null(control[["lambda0",exact=TRUE]])) { control$lambda0 <- NULL }
  #Which method to use for computing eq3a
  if (is.null(control[["eq3a.method",exact=TRUE]])) { control$eq3a.method <- "R" } else {
    control$eq3a.method <- match.arg(control$eq3a.method,c("R","C"))
  }
  #Hook function definition
  if (is.null(control[["hookFun",exact=TRUE]])) { control$hookFun <- function(Y,lambda,...) {} }

  #If the eps and iter.max arguments are too short, make them length 2.
  if (length(control$eps)==1) control$eps <- rep(control$eps,2)
  if (length(control$iter.max)==1) control$iter.max <- rep(control$iter.max,2)
  
  #Compute the estimate to report (i.e. use 2nd component of the args)
  if (control$verbose) {
    cat("Back-projecting with k=",control$k," to get lambda estimate.\n")
  }
  stsk <- backprojNP.fit(sts, incu.pmf=incu.pmf,k=control$k,eps=control$eps[2],iter.max=control$iter.max[2],verbose=control$verbose,lambda0=control$lambda0,hookFun=control$hookFun,eq3a.method=control$eq3a.method)
  #Fix control slot
  stsk@control <- control 

  #If no bootstrap to do return object right away as stsBP object.
  if (control$B<=0) {
    if (control$verbose) { cat("No bootstrap CIs calculated as requested.\n") }
    stsk <- as(stsk,"stsBP")
    return(stsk)
  }

  #Call back-project function without smoothing, i.e. with k=0.
  if (control$verbose) {
    cat("Back-projecting with k=",0," to get lambda estimate for parametric bootstrap.\n")
  }
  sts0 <- backprojNP.fit(sts, incu.pmf=incu.pmf,k=0,eps=control$eps[1],iter.max=control$iter.max[1],verbose=control$verbose,lambda0=control$lambda0,hookFun=control$hookFun, eq3a.method=control$eq3a.method)

  ###########################################################################
  #Create bootstrap samples and loop for each sample while storing the result
  ###########################################################################
  sts.boot <- sts0
  #Define object to return
  lambda <- array(NA,dim=c(nrow(sts),ncol(sts),control$B))

  #Define PMF of incubation time which does safe handling of values
  #outside the support of the incubation time.
  dincu <- function(x,i) {
    notInSupport <- x<0 | x>=length(incu.pmf[,i])
    #Give index -1 to invalid queries
    x[notInSupport] <- -1
    return(c(0,incu.pmf[,i])[x+2])
  }
  
  #Loop in order to create the sample
  for (b in 1:control$B) {
    if (control$verbose) { cat("Bootstrap sample ",b,"/",control$B,"\n") }
    
    #Compute convolution for the mean of the observations 
    mu <- matrix(0, nrow=nrow(sts0), ncol=ncol(sts0))
    #Perform the convolution for each series
    for (i in 1:ncol(sts)) {
      for (t in 1:nrow(mu)) {
        for (s in 0:(t-1)) {
          mu[t,i] <- mu[t,i] + upperbound(sts0)[t-s,i] * dincu(s,i)
        }
      }
    }
    
    #Create new observations in the observed slot.
    observed(sts.boot) <- matrix(rpois(prod(dim(sts.boot)),lambda=mu),ncol=ncol(sts0))

    #Run the backprojection on the bootstrap sample. Use original result
    #as starting value.
    sts.boot <- backprojNP.fit(sts.boot, incu.pmf=incu.pmf,k=control$k,eps=control$eps[2],iter.max=control$iter.max[2],verbose=control$verbose,lambda0=upperbound(stsk),hookFun=control$hookFun, eq3a.method=control$eq3a.method)
    #Extract the result of the b'th backprojection
    lambda[,,b] <- upperbound(sts.boot)
  }

  #Compute an equal tailed (1-alpha)*100% confidence intervals based on the
  #bootstrap samples. The dimension is (ci.low,ci.high) x time x series
  ci <- apply(lambda,MARGIN=c(1,2), quantile, p=c(control$alpha/2,1-control$alpha/2)) 

  #Convert output to stsBP object and add information to the extra slots
  stsk <- as(stsk,"stsBP")

  #Add extra slots
  stsk@ci <- ci
  stsk@lambda <- lambda
  stsk@control <- control

  #Done
  return(stsk)
}

