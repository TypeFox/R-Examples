#' @include generics.R
#' @include Class-Weight.R
#' @include Class-LagOperator.R
#' 
NULL
################################################################################
#' Class for a lag-window type estimator.
#'
#' For a given time series Y a lag-window estimator of the Form
#' \deqn{\hat{f}(\omega) = \sum_{|k|< n-1 } K_n(k) \Gamma(Y_0,Y_k) \exp(-i \omega k)} 
#' will be calculated on initalization. The \code{LagKernelWeight} K_n is determined
#' by the slot \code{weight} and the \code{LagOperator} \eqn{\Gamma(Y_0,Y_k)} is defined 
#' by the slot lagOp.
#' 
#' Currently, the implementation of this class allows only for the analysis of
#' univariate time series. 
#'
#' @name LagEstimator-class
#' @aliases LagEstimator
#' @exportClass LagEstimator
#'
#' @keywords S4-classes
#'
#' @slot Y the time series where the lag estimator was applied one
#' @slot weight a \code{\link{Weight}} object to be used as lag window
#' @slot lagOp a \code{\link{LagOperator}} object that determines which
#'             kind of bivariate structure should be calculated.
#' @slot env An environment to allow for slots which need to be
#'       accessable in a call-by-reference manner:
#'       \describe{
#'         \item{\code{sdNaive}}{An array used for storage of the naively
#'             estimated standard deviations of the smoothed periodogram.}
#'         \item{\code{sdNaive.done}}{a flag indicating whether \code{sdNaive}
#'             has been set yet.}}
#'
################################################################################
setClass(
  Class = "LagEstimator",
  representation=representation(
    env = "environment",
    Y = "numeric",
    weight = "Weight",
    lagOp = "LagOperator"
  ),
  contains = "QSpecQuantity"
)

#' @importFrom stats quantile
#' @importFrom stats mvfft
setMethod(
  f = "initialize",
  signature = "LagEstimator",
  definition = function(.Object, Y, frequencies, lagOp, weight){
    
     levels.1 = lagOp@levels.1
     levels.2 = lagOp@levels.2
    
    .Object@Y = Y
    .Object@frequencies = frequenciesValidator(frequencies, length(Y), steps=1:6)
    .Object@levels[[1]] = levels.1
    .Object@levels[[2]] = levels.2
    .Object@weight = weight
    .Object@lagOp = lagOp
    .Object@env <- new.env(parent=emptyenv())

    .Object@env$sdNaive.done <- FALSE
    .Object@env$sdNaive <- array()
    .Object@env$sdBoot.done <- FALSE
    .Object@env$sdBoot <- array()
    
    levels.all = union(levels.1, levels.2)
    ln1 = length(levels.1)
    ln2 = length(levels.2)
    ln = length(levels.all)
    Q <- Y
    Q <- quantile(Y, probs = levels.all)
    
    B <- getB(lagOp)
    
    res <- array(dim = c(length(Y), ln1, ln2, B+1))
    
    Fhat = getValues(lagOp, levels.1 = levels.1, levels.2 = levels.2)
    n = dim(Fhat)[1]
    Kernel = getValues(weight)
    if(length(Kernel)<n){
      warning("number of 'weights' is too small, filled up with zeroes")
      Kernel = c(Kernel,rep(0, n - length(Kernel)))
    }
    if(length(Kernel) > n){
      Kernel = Kernel[1:n]
    }
    
    for (l1 in 1:ln1) {
      for (l2 in 1:ln2) {  #\tau_2,\tau_1 can be calculated from \tau_1,\tau_2 
        if((l1 > l2)&&!(levels.1[l1]==levels.2[l2])&&is.element(levels.1[l1],levels.2)&&is.element(levels.2[l2],levels.1)){ 
          res[,l1,l2,] <- Conj(res[,l2,l1,])
          next
        }
        A1 <- matrix(Fhat[,l1,l2,], ncol = B+1)
        A1 <- A1# - levels.1[l1] * levels.2[l2]
        A1 <- A1*Kernel
        
        A2 <- matrix(Fhat[,l2,l1,], ncol = B+1)
        A2 <- A2# - levels.1[l2] * levels.2[l1]  				 # same as above!
        A2 <- A2*Kernel
        
        A0 <- matrix(rep(Fhat[1,l2,l1,], nrow(A1)), ncol = B+1, byrow = TRUE)
        A0 <- A0# - levels.1[l2] * levels.2[l1]					 # same as above!
        A0 <- A0*Kernel[1]
        res[,l1,l2,] <- 1/(2*pi) * (mvfft(A1) + Conj(mvfft(A2)) - A0)
      }
    }
    .Object@values = res
    return(.Object)
  })

################################################################################
#' Get values from a lag-window type estimator.
#'
#' The returned array of \code{values} is of dimension \code{[J,K1,K2]},
#' where \code{J=length(frequencies)}, \code{K1=length(levels.1)} and
#' \code{K2=length(levels.2))}. At position \code{(j,k1,k2)}
#' the returned value is the one corresponding to \code{frequencies[j]},
#' \code{levels.1[k1]} and \code{levels.2[k2]} that are closest to the
#' \code{frequencies}, \code{levels.1} and \code{levels.2}
#' available in \code{object}; \code{\link{closest.pos}} is used to determine
#' what closest to means. 
#' 
#' @name getValues-LagEstimator
#' @aliases getValues,LagEstimator-method
#'
#' @keywords Access-functions
#'
#' @param object \code{LagEstimator} of which to get the values
#' @param frequencies a vector of frequencies for which to get the values
#' @param levels.1 the first vector of levels for which to get the values
#' @param levels.2 the second vector of levels for which to get the values
#'
#' @return Returns data from the array \code{values} that's a slot of
#'          \code{object}.
#'
#' @seealso
#' An example on how to use this function is analogously to the example given in
#' \code{\link{getValues-QuantilePG}}.
################################################################################

setMethod(f = "getValues",
          signature = signature(
            "LagEstimator"),
          definition = function(object,
                                frequencies=2*pi*(0:(length(object@Y)-1))/length(object@Y),
                                levels.1=getLevels(object,1),
                                levels.2=getLevels(object,2)) {
            
            # workaround: default values don't seem to work for generic functions?
            if (!hasArg(frequencies)) {
              frequencies <- 2*pi*(0:(length(object@Y)-1))/length(object@Y)
            }
            if (!hasArg(levels.1)) {
              levels.1 <- object@levels[[1]]
            }
            if (!hasArg(levels.2)) {
              levels.2 <- object@levels[[2]]
            }
            # end: workaround
            
            ##############################
            ## (Similar) Code also in Class-FreqRep!!!
            ## (Similar) Code also in Class-QuantileSD!!!
            ##############################
            
            # Transform all frequencies to [0,2pi)
            frequencies <- frequencies %% (2*pi)
            
            # Create an aux vector with all available frequencies
            oF <- object@frequencies
            f <- frequencies
            
            # returns TRUE if x c y
            subsetequal.approx <- function(x,y) {
              X <- round(x, .Machine$double.exponent-2)
              Y <- round(y, .Machine$double.exponent-2)
              return(setequal(X,intersect(X,Y)))
            }
            
            C1 <- subsetequal.approx(f[f <= pi], oF)
            C2 <- subsetequal.approx(f[f > pi], 2*pi - oF[which(oF != 0 & oF != pi)])
            
            if (!(C1 & C2)) {
              warning("Not all 'values' for 'frequencies' requested were available. 'values' for the next available Fourier frequencies are returned.")
            }
            
            # Select columns
            c.1.pos <- closest.pos(object@levels[[1]],levels.1)
            c.2.pos <- closest.pos(object@levels[[2]],levels.2)
            
            if (!subsetequal.approx(levels.1, object@levels[[1]])) {
              warning("Not all 'values' for 'levels.1' requested were available. 'values' for the next available level are returned.")
            }
            
            if (!subsetequal.approx(levels.2, object@levels[[2]])) {
              warning("Not all 'values' for 'levels.2' requested were available. 'values' for the next available level are returned.")
            }
            
            
            J <- length(frequencies)
            K1 <- length(levels.1)
            K2 <- length(levels.2)
            res <- array(dim=c(J, K1, K2, getB(object@lagOp)+1))
            
            
            if (class(object@weight) == "LagKernelWeight") {
              
              # Select rows
              r1.pos <- closest.pos(oF, f[f <= pi])
              r2.pos <- closest.pos(-1*(2*pi-oF),-1*f[f > pi])
              
              if (length(r1.pos) > 0) {
                res[which(f <= pi),,,] <- object@values[r1.pos,c.1.pos,c.2.pos,]
              }
              if (length(r2.pos) > 0) {
                res[which(f > pi),,,] <- Conj(object@values[r2.pos,c.1.pos,c.2.pos,])
              }
              
            }
            
            return(res)
          }
)

################################################################################
#' Get estimates for the standard deviation of the lagEstimator derived from
#' the asymptotics (see Birr et al (2015))
#'
#' Determines and returns an array of dimension \code{[J,K1,K2]},
#' where \code{J=length(frequencies)}, \code{K1=length(levels.1)}, and
#' \code{K2=length(levels.2))}.
#' At position \code{(j,k1,k2)} the returned value is the standard deviation 
#' estimated corresponding to \code{frequencies[j]}, \code{levels.1[k1]} and 
#' \code{levels.2[k2]} that are closest to the \code{frequencies}, \code{levels.1}
#' and \code{levels.2} available in \code{object}; \code{\link{closest.pos}} is 
#' used to determine what closest to means.
#'
#' Requires that the \code{\link{LagEstimator}} is available at all Fourier
#' frequencies from \eqn{(0,\pi]}{(0,pi]}. If this is not the case the missing
#' values are imputed by taking one that is available and has a frequency
#' that is closest to the missing Fourier frequency; \code{closest.pos} is used
#' to determine which one this is.
#'
#' Note the ``standard deviation'' estimated here is not the square root of the
#' complex-valued variance. It's real part is the square root of the variance
#' of the real part of the estimator and the imaginary part is the square root
#' of the imaginary part of the variance of the estimator.
#'
#' @name getSdNaive-LagEstimator
#' @aliases getSdNaive,LagEstimator-method
#' 
#' @importFrom stats integrate
#'
#' @keywords Access-functions
#'
#' @param object \code{\link{LagEstimator}} of which to get the estimates for the
#'                standard deviation.
#' @param frequencies a vector of frequencies for which to get the result
#' @param levels.1 the first vector of levels for which to get the result
#' @param levels.2 the second vector of levels for which to get the result
#' 
#'
#' @return Returns the estimate described above.
#'
################################################################################
setMethod(f = "getSdNaive",
          signature = signature(
            object = "LagEstimator"),
          definition = function(object,
                                frequencies=2*pi*(0:(length(object@Y)-1))/length(object@Y),
                                levels.1=getLevels(object,1),
                                levels.2=getLevels(object,2)) {
            
            if (class(object@weight) != "LagKernelWeight") {
              stop("getSdNaive currently only available for 'LagKernelWeight'.")
            }
            
            
            objLevels1 <- object@levels[[1]]
            objLevels2 <- object@levels[[2]]
            Y <- object@Y
            # workaround: default values don't seem to work for generic functions?
            if (!hasArg(frequencies)) {
              frequencies <- 2*pi*(0:(length(Y)-1))/length(Y)
            }
            if (!hasArg(levels.1)) {
              levels.1 <- objLevels1
            }
            if (!hasArg(levels.2)) {
              levels.2 <- objLevels2
            }
            # end: workaround
            
            if (object@env$sdNaive.done == FALSE) {
              
              weight <- object@weight
              Bn = weight@bw
              N <- length(Y)
              K1 <- length(objLevels1)
              K2 <- length(objLevels2)
               
                V <- getValues(object, frequencies = 2*pi*(1:(N-1))/N)    
                Sigma <- array(0,dim=c(N-1,K1,K2))
                K_Int = integrate(function(x){getW(weight)(x)^2},lower = -1, upper = 1)$value
                #print(K_Int)
                for(k1 in 1:K1){
                  for(k2 in k1:K2){
                    Sigma[,k1,k2] = sqrt(K_Int*Re(V[,k1,k1,1])*Re(V[,k2,k2,1]))
                    Sigma[,k2,k1] = Sigma[,k1,k2]
                  }
                }
                
              object@env$sdNaive.done <- TRUE
              object@env$sdNaive <- sqrt(Bn/N)*(Re(Sigma)>0)*(complex(real = Re(Sigma),imaginary = Re(Sigma)))
            } # End of 'if (object@env$sdNaive.done == FALSE) {'
            
            resObj <- object@env$sdNaive
            ##############################
            ## (Similar) Code also in Class-FreqRep!!!
            ##############################
            
            # Transform all frequencies to [0,2pi)
            frequencies <- frequencies %% (2*pi)
            
            # Create an aux vector with all available frequencies
            oF <- object@frequencies
            f <- frequencies
            
            # returns TRUE if x c y
            subsetequal.approx <- function(x,y) {
              X <- round(x, .Machine$double.exponent-2)
              Y <- round(y, .Machine$double.exponent-2)
              return(setequal(X,intersect(X,Y)))
            }
            
            C1 <- subsetequal.approx(f[f <= pi], oF)
            C2 <- subsetequal.approx(f[f > pi], 2*pi - oF[which(oF != 0 & oF != pi)])
            
            # Ist dann der Fall wenn die frequencies nicht geeignete FF sind!!
            if (!(C1 & C2)) {
              warning("Not all 'values' for 'frequencies' requested were available. 'values' for the next available Fourier frequencies are returned.")
            }
            
            # Select columns
            c.1.pos <- closest.pos(objLevels1,levels.1)
            c.2.pos <- closest.pos(objLevels2,levels.2)
            
            # Select rows
            r1.pos <- closest.pos(oF,f[f <= pi])
            r2.pos <- closest.pos(-1*(2*pi-oF),-1*f[f > pi])
            
            J <- length(frequencies)
            K1 <- length(levels.1)
            K2 <- length(levels.2)
            res <- array(dim=c(J, K1, K2))
            
            if (length(r1.pos) > 0) {
              res[which(f <= pi),,] <- resObj[r1.pos,c.1.pos,c.2.pos]
            }
            if (length(r2.pos) > 0) {
              res[which(f > pi),,] <- Conj(resObj[r2.pos,c.1.pos,c.2.pos])
            }
            
            return(res)
          }
)

################################################################################
#' Get bootstrap estimates for the standard deviation of the lag-window type
#' estimator.
#'
#' Determines and returns an array of dimension \code{[J,K1,K2]},
#' where \code{J=length(frequencies)}, \code{K1=length(levels.1)}, and
#' \code{K2=length(levels.2))}.
#' At position \code{(j,k1,k2)} the real part of the returned value is the
#' standard deviation estimated from the real parts of the bootstrap
#' replications and the imaginary part of the returned value is the standard
#' deviation estimated from the imaginary part of the bootstrap replications.
#' The estimate is determined from those bootstrap replicates of the estimator
#' that have
#' \code{frequencies[j]}, \code{levels.1[k1]} and \code{levels.2[k2]} closest
#' to the \code{frequencies}, \code{levels.1} and \code{levels.2}
#' available in \code{object}; \code{\link{closest.pos}} is used to determine
#' what closest to means.
#'
#' Requires that the \code{\link{LagEstimator}} is available at all Fourier
#' frequencies from \eqn{(0,\pi]}{(0,pi]}. If this is not the case the missing
#' values are imputed by taking one that is available and has a frequency
#' that is closest to the missing Fourier frequency; \code{closest.pos} is used
#' to determine which one this is.
#'
#' If there are no bootstrap replicates available (i. e., \code{B == 0}) an
#' error is returned.
#'
#' Note the ``standard deviation'' estimated here is not the square root of the
#' complex-valued variance. It's real part is the square root of the variance
#' of the real part of the estimator and the imaginary part is the square root
#' of the imaginary part of the variance of the estimator.
#'
#' @name getSdBoot-LagEstimator
#' @aliases getSdBoot,LagEstimator-method
#' 
#' @importFrom stats sd
#'
#' @keywords Access-functions
#'
#' @param object \code{\link{LagEstimator}} of which to get the bootstrap
#' 								estimates for the standard deviation.
#' @param frequencies a vector of frequencies for which to get the result
#' @param levels.1 the first vector of levels for which to get the result
#' @param levels.2 the second vector of levels for which to get the result
#'
#' @return Returns the estimate described above.
################################################################################

setMethod(f = "getSdBoot",
    signature = signature(
        object = "LagEstimator"),
    definition = function(object,
        frequencies=2*pi*(0:(length(object@lagOp@Y)-1))/length(object@lagOp@Y),
        levels.1=getLevels(object,1),
        levels.2=getLevels(object,2)) {
      
      if (class(getWeight(object)) != "LagKernelWeight") {
        stop("getSdBoot currently only available for 'LagKernelWeight'.")
      }
      
      # workaround: default values don't seem to work for generic functions?
      if (!hasArg(frequencies)) {
        frequencies <- 2*pi*(0:(length(object@lagOp@Y)-1))/length(object@lagOp@Y)
      }
      if (!hasArg(levels.1)) {
        levels.1 <- object@levels[[1]]
      }
      if (!hasArg(levels.2)) {
        levels.2 <- object@levels[[2]]
      }
      # end: workaround
      
      #if (object@env$sdBoot.done) {
        
        complex.var <- function(x) {
          return(complex(real = sd(Re(x)), imaginary = sd(Im(x))))
        }
        
        #N <- length(object@qPG@freqRep@Y)
        #K1 <- length(object@levels[[1]])
        #K2 <- length(object@levels[[2]])
        B <- getB(object@lagOp)
        
        v <- getValues(object, frequencies = frequencies,
            levels.1 = levels.1, levels.2 = levels.2)[,,,2:(B+1), drop=F]
        object@env$sdBoot <- apply(v, c(1,2,3), complex.var)
      #  object@env$sdBoot.done <- TRUE
      #}
      return(object@env$sdBoot)
    }
)

################################################################################
#' Get pointwise confidence intervals for the quantile spectral density kernel
#'
#' Returns a list of two arrays \code{lowerCIs} and \code{upperCIs} that contain
#' the upper and lower limits for a level \code{1-alpha} confidence interval of
#' the copula spectral density kernel. Each array is of dimension \code{[J,K1,K2]},
#' where \code{J=length(frequencies)}, \code{K1=length(levels.1)}, and
#' \code{K2=length(levels.2))}.
#' At position \code{(j,k1,k2)} the real (imaginary) part of the returned values
#' are the bounds of the confidence interval for the the real (imaginary) part
#' of the quantile spectrum, which corresponds to
#' \code{frequencies[j]}, \code{levels.1[k1]} and \code{levels.2[k2]} closest
#' to the Fourier frequencies, \code{levels.1} and \code{levels.2}
#' available in \code{object}; \code{\link{closest.pos}} is used to determine
#' what closest to means.
#'
#' Currently, only one \code{type} of confidence interval is
#' available:
#' \itemize{
#'   \item \code{"naive.sd"}: confidence intervals based on the asymptotic
#'           normality of the lag-window estimator; standard deviations
#'           are estimated using \code{\link{getSdNaive}}.
#'  }
#'
#' @name getPointwiseCIs-LagEstimator
#' @aliases getPointwiseCIs,LagEstimator-method
#'
#' @importFrom stats qnorm
#' @importFrom stats quantile
#' 
#' @keywords Access-functions
#'
#' @param object \code{LagEstimator} of which to get the confidence intervals
#' @param frequencies a vector of frequencies for which to get the result
#' @param levels.1 the first vector of levels for which to get the result
#' @param levels.2 the second vector of levels for which to get the result
#' @param alpha the level of the confidence interval; must be from \eqn{(0,1)}
#' @param type a flag indicating which type of confidence interval should be
#'         returned; can only take one values at the moment.
#'
#' @return Returns a named list of two arrays \code{lowerCIS} and \code{upperCIs}
#'          containing the lower and upper bounds for the confidence intervals.
#'
#' @examples
#' lagEst <- lagEstimator(rnorm(2^10), levels.1=0.5)
#' CI.upper <- Re(getPointwiseCIs(lagEst)$upperCIs[,1,1])
#' CI.lower <- Re(getPointwiseCIs(lagEst)$lowerCIs[,1,1])
#' freq = 2*pi*(0:1023)/1024
#' plot(x = freq, y = rep(0.25/(2*pi),1024),
#'    ylim=c(min(CI.lower), max(CI.upper)),
#'    type="l", col="red") # true spectrum
#' lines(x = freq, y = CI.upper)
#' lines(x = freq, y = CI.lower)
################################################################################
setMethod(f = "getPointwiseCIs",
    signature = signature(
      object = "LagEstimator"),
    definition = function(object,
                          frequencies=2*pi*(0:(length(object@Y)-1))/length(object@Y),
                          levels.1=getLevels(object,1),
                          levels.2=getLevels(object,2),
                          alpha=.1, type=c("naive.sd", "boot.sd", "boot.full")) {
      
      # workaround: default values don't seem to work for generic functions?
      if (!hasArg(frequencies)) {
        frequencies <- 2*pi*(0:(length(object@Y)-1))/length(object@Y)
      }
      if (!hasArg(levels.1)) {
        levels.1 <- object@levels[[1]]
      }
      if (!hasArg(levels.2)) {
        levels.2 <- object@levels[[2]]
      }
      if (!hasArg(alpha)) {
        alpha <- 0.1
      }
      if (!hasArg(type)) {
        type <- "naive.sd"
      }
      # end: workaround
      
      type <- match.arg(type)[1]
      switch(type,
             "naive.sd" = {
               sdEstim <- getSdNaive(object,
                                     frequencies = frequencies,
                                     levels.1 = levels.1,
                                     levels.2 = levels.2)},
             "boot.sd" = {
               sdEstim <- getSdBoot(object,
                                    frequencies = frequencies,
                                    levels.1 = levels.1,
                                    levels.2 = levels.2)}
      )
      
      J <- length(frequencies)
      K1 <- length(levels.1)
      K2 <- length(levels.2)
      
      
      if (type == "naive.sd" || type == "boot.sd") {
        v <- getValues(object,
                       frequencies = frequencies,
                       levels.1 = levels.1,
                       levels.2 = levels.2)[,,,1]
        upperCIs <- array(v + sdEstim * qnorm(1-alpha/2), dim = c(J, K1, K2))
        lowerCIs <- array(v + sdEstim * qnorm(alpha/2), dim = c(J, K1, K2))
      } else if (type == "boot.full") {
        # TODO: Error Msg ausgeben falls B == 0
        B <- getB(object@lagOp)
        v <- getValues(object,
           frequencies = frequencies,
           levels.1 = levels.1,
           levels.2 = levels.2)[,,,2:(B+1), drop=FALSE]
        uQuantile <- function(x) {complex(real = quantile(Re(x),1-alpha/2),
              imaginary = quantile(Im(x),1-alpha/2))}
        lQuantile <- function(x) {complex(real = quantile(Re(x),alpha/2),
              imaginary = quantile(Im(x),alpha/2))}
        
        upperCIs <- apply(v, c(1,2,3), uQuantile)
        lowerCIs <- apply(v, c(1,2,3), lQuantile)
      }
      
      res <- list(lowerCIs = lowerCIs, upperCIs = upperCIs)
      return(res)
    }
)


  ################################################################################
  #' Create an instance of the \code{LagEstimator} class.
  #'
  #' A \code{LagEstimator} object can be created from \code{numeric}, a \code{ts},
  #' or a \code{zoo} object. Also a \code{\link{LagOperator}} and a 
  #' \code{\link{Weight}} object can be used to create different types of 
  #' estimators.
  #'
  #'
  #' @name LagEstimator-constructor
  #' @aliases lagEstimator
  #' @export
  #'
  #' @keywords Constructors
  #'
  #' @param Y a time series (\code{numeric}, \code{ts}, or \code{zoo} object) or a 
  #'          \code{\link{LagOperator}} from which to determine the \code{LagEstimator}
  #' @param frequencies A vector containing (Fourier-)frequencies at which to determine the
  #'                     smoothed periodogram.
  #' @param levels.1 the first vector of levels for which to compute the LagEstimator
  #' @param levels.2 the second vector of levels for which to compute the LagEstimator               
  #' @param weight Object of type \code{\link{Weight}} to be used for smoothing.
  #' @param type if \code{Y} is a time series, this indicates which LagOperator will be used 
  #' @return Returns an instance of \code{LagEstimator}.
  #'
  #' @examples
  #' Y <- rnorm(100)
  #' levels.1 <- c(0.1,0.5,0.9)
  #' weight <- lagKernelWeight(W = WParzen,  bw = 10, K = length(Y))
  #' lagOp <- clippedCov(Y,levels.1 = levels.1)
  #' lagEst <- lagEstimator(lagOp, weight = weight)
################################################################################
lagEstimator <- function(Y,
                         frequencies=2*pi/length(Y) * 0:(length(Y)-1),
                         levels.1 = .5,
                         levels.2 = levels.1,
                         weight = lagKernelWeight(K = length(Y),bw = 100),
                         type = c("clippedCov")
                         ){
  
  #ToDo Check if params are okay
  #     Transform to Fourier frequencies
  if (class(Y) == "numeric") {
    versConstr <- 1
    Y <- Y
  } else if (class(Y) == "ts") {
    versConstr <- 1
    Y <- Y[1:length(Y)]
  } else if (class(Y) == "zoo") {
    versConstr <- 1
    Y <- coredata(Y)
  } else if ( is(Y,"LagOperator")) {
    lagOp = Y
    versConstr <- 2
    if (!hasArg(frequencies)) {
      Z <- Y@Y
      frequencies <- 2*pi/length(Z) * 0:(length(Z)-1)
    }
    
    if (!hasArg(levels.1)) {
      levels.1 <- getLevels(Y,1)
    }
    if (!hasArg(levels.2)) {
      levels.2 <- getLevels(Y,2)
    }
    Y = Z
  } else {
    stop("object is neither 'numeric', 'ts', 'zoo', nor 'LagOperator'.")
  }
  if(versConstr == 1){
    if (type == "clippedCov"){
      lagOp = clippedCov(Y, levels.1 = levels.1, levels.2 = levels.2)  
    }
    
  }

   return(new(Class = "LagEstimator",
           Y = Y,
           frequencies = frequencies,
           weight = weight,
           lagOp = lagOp)
   )
  
  
}

################################################################################
#' Get associated \code{\link{Weight}} from a \code{\link{LagEstimator}}.
#'
#' @name getWeight-LagEstimator
#' @aliases getWeight,LagEstimator-method
#'
#' @keywords Access-association-functions
#'
#' @param object \code{LagEstimator} from which to get the \code{Weight}.
#' @return Returns the \code{\link{Weight}} object associated.
################################################################################
setMethod(f = "getWeight",
    signature = "LagEstimator",
    definition = function(object) {
      return(object@weight)
    }
)

################################################################################
#' Get associated \code{\link{LagOperator}} from a \code{\link{LagEstimator}}.
#'
#' @name getLagOperator-LagEstimator
#' @aliases getLagOperator,LagEstimator-method
#'
#' @keywords Access-association-functions
#'
#' @param object \code{LagEstimator} from which to get the \code{LagOperator}.
#' @return Returns the \code{\link{LagOperator}} object associated.
################################################################################
setMethod(f = "getLagOperator",
    signature = "LagEstimator",
    definition = function(object) {
      return(object@lagOp)
    }
)
  
################################################################################
#' Plot the values of a \code{\link{LagEstimator}}.
#'
#' Creates a \code{K} x \code{K} plot displaying all levels combinations from the
#' argument \code{levels}.  
#' In each of the subplots either the real part (on and below the diagonal;
#' i. e., \eqn{\tau_1 \leq \tau_2}{tau1 <= tau2}) or the imaginary parts
#' (above the diagonal; i. e., \eqn{\tau_1 > \tau_2}{tau1 > tau2}) of
#' the lag-window estimator, for the combination of levels \eqn{\tau_1}{tau1}
#'  and \eqn{\tau_2}{tau2} denoted on the left and bottom margin of the plot are displayed.
#'
#' @name plot-LagEstimator
#' @aliases plot,LagEstimator,ANY-method
#' @export
#'
#' @importFrom abind abind
#'
#' @param x  The \code{\link{LagEstimator}} object to plot
#' @param ptw.CIs the confidence level for the confidence intervals to be
#'                 displayed; must be a number from [0,1]; if null, then no
#'                 confidence intervals will be plotted.
#'                 
#' @param ratio quotient of width over height of the subplots; use this
#'               parameter to produce landscape or portrait shaped plots.
#' @param widthlab width for the labels (left and bottom); default is
#'                  \code{lcm(1)}, cf. \code{\link[graphics]{layout}}.
#' @param xlab label that will be shown on the bottom of the plots; can be
#'              an expression (for formulas), characters or \code{NULL} to
#'              force omission (to save space).
#' @param ylab label that will be shown on the left side of the plots;
#'              can be an expression (for formulas), characters or
#'              \code{NULL} to force omission (to save space).
#' @param type.scaling a method for scaling of the subplots; currently there
#'                      are three options: \code{"individual"} will scale each of the
#'                      \code{K^2} subplots to minimum and maximum of the values
#'                      in that plot, \code{"real-imaginary"} will scale each of the
#'                      subplots displaying real parts and each of the subplots
#'                      displaying imaginary parts to the minimum and maximum of
#'                      the values display in these subportion of plots. The
#'                      option \code{"all"} will scale the subplots to the minimum and
#'                      maximum in all of the subplots.
#' @param type.CIs indicates the method to be used for determining the
#'                 confidence intervals; the methods available are those
#'                  provided by
#'                 \code{\link{getPointwiseCIs-LagEstimator}}.
#' @param frequencies a set of frequencies for which the values are to be
#'                    plotted.
#' @param levels a set of levels for which the values are to be plotted.
#'
#' @return Returns the plot described in the Description section.
#' 
#' See Birr et al. (2015)
#' @references
#' Birr, S., Volgushev, S., Kley, T., Dette, H. & Hallin, M. (2015).
#' Quantile Spectral Analysis for Locally Stationary Time Series.
#' \url{http://arxiv.org/abs/1404.4605}.
################################################################################

setMethod(f = "plot",
          signature = signature("LagEstimator"),
          definition = function(x,
                                ptw.CIs = 0.1, 
                                ratio = 3/2, widthlab = lcm(1), xlab = expression(omega/2*pi), ylab = NULL,
                                type.scaling = c("individual", "real-imaginary", "all"),
                                frequencies=x@frequencies,
                                type.CIs = c("naive.sd"),
                                levels=intersect(x@levels[[1]], x@levels[[2]])) {
            
            def.par <- par(no.readonly = TRUE) # save default, for resetting...
            
            # workaround: default values don't seem to work for generic functions?
            if (!hasArg(ptw.CIs)) {
              ptw.CIs <- 0.1
            }
            if (!hasArg(frequencies)) {
              frequencies <- x@frequencies
            }
            if (!hasArg(levels)) {
              levels <- intersect(x@levels[[1]], x@levels[[2]])
            }
            if (!hasArg(ptw.CIs)) {
              ptw.CIs <- 0
            }
            if (!hasArg(ratio)) {
              ratio <- 3/2
            }
            if (!hasArg(widthlab)) {
              widthlab <- lcm(1)
            }
            if (!hasArg(xlab)) {
              xlab <- expression(omega/2*pi)
            }
            if (!hasArg(ylab)) {
              ylab <- NULL
            }
            if (!hasArg(type.scaling)) {
              type.scaling <- c("individual", "real-imaginary", "all")
            }
            # end: workaround
            
            if (length(levels) == 0) {
              stop("There has to be at least one level to plot.")
            }
            
            tryCatch({
              
              K <- length(levels)
              values <- getValues(x, frequencies = frequencies,
                                  levels.1=levels, levels.2=levels)
            #text.headline <- x@weight@descr
              if (ptw.CIs > 0) {
                CI <- getPointwiseCIs(x, frequencies = frequencies,
                                      alpha=ptw.CIs, type=type.CIs,
                                      levels.1=levels, levels.2=levels)
                lowerCIs  <- CI$lowerCIs
                upperCIs  <- CI$upperCIs
                #text.headline <- (paste(text.headline, ", includes ",1-ptw.CIs,"-CI (ptw. of type '",type.CIs,"')",sep=""))
              }
              
              X <- frequencies/(2*pi)
              
              allVals <- array(values[,,,1], dim=c(length(X), K, K))
             
              if (ptw.CIs > 0) {
                allVals <- abind(allVals, lowerCIs, upperCIs, along=1)
              }
              type.scaling <- match.arg(type.scaling)[1]
              
              p <- K
              M <- matrix(1:p^2, ncol=p)
              M.heights <- rep(1,p)
              M.widths  <- rep(ratio,p)
              
              # Add places for tau labels
              M <- cbind((p^2+1):(p^2+p),M)
              M.widths <- c(widthlab,M.widths)
              M <- rbind(M,c(0,(p^2+p+1):(p^2+2*p)))
              M.heights <- c(M.heights, widthlab)
              
              i <- (p^2+2*p+1)
              # Add places for labels
              if (length(xlab)>0) {
                M.heights <- c(M.heights, widthlab)
                M <- rbind(M,c(rep(0,length(M.widths)-p),rep(i,p)))
                i <- i + 1
              }
              
              if (length(ylab)>0) {
                M <- cbind(c(rep(i,p),rep(0,length(M.heights)-p)),M)
                M.widths <- c(widthlab,M.widths)
              }
              
              nf <- layout(M, M.widths, M.heights, TRUE)
              
              par(mar=c(2,2,1,1))
              
              for (i1 in 1:K) {
                for (i2 in 1:K) {
                  if (i2 >= i1) {
                    switch(type.scaling,
                           "individual" = {
                             y.min <- min(Re(allVals[,i1,i2]))
                             y.max <- max(Re(allVals[,i1,i2]))},
                           "real-imaginary" = {
                             y.min <- min(Re(allVals))
                             y.max <- max(Re(allVals))},
                           "all" = {
                             y.min <- min(Re(allVals),Im(allVals))
                             y.max <- max(Re(allVals),Im(allVals))}
                    )
                    
                    plot(x=0,y=0, type="n", xlab="", ylab="", #xlab=xl, ylab=yl,
                         xlim=c(min(X), max(X)), ylim=c(y.min, y.max))
                    if (ptw.CIs > 0) {
                      polygon(x=c(X,rev(X)), y=c(Re(lowerCIs[,i1,i2]),rev(Re(upperCIs[,i1,i2]))),
                              col="lightgray", border=NA)
                    }
                    lines(x=X, y=Re(values[,i1,i2,1]),
                          ylim=c(min(Re(allVals)), max(Re(allVals))),
                          type="l", col="blue")
                  } else {
                    switch(type.scaling,
                           "individual" = {
                             y.min <- min(Im(allVals[,i1,i2]))
                             y.max <- max(Im(allVals[,i1,i2]))},
                           "real-imaginary" = {
                             y.min <- min(Im(allVals))
                             y.max <- max(Im(allVals))},
                           "all" = {
                             y.min <- min(Re(allVals),Im(allVals))
                             y.max <- max(Re(allVals),Im(allVals))}
                    )
                    plot(x=0,y=0, type="n", xlab="", ylab="", #xlab=xl, ylab=yl,
                         xlim=c(min(X), max(X)), ylim=c(y.min, y.max))
                    if (ptw.CIs > 0) {
                      polygon(x=c(X,rev(X)), y=c(Im(lowerCIs[,i1,i2]),rev(Im(upperCIs[,i1,i2]))),
                              col="lightgray", border=NA)
                    }
                    lines(x=X, y=Im(values[,i1,i2,1]),
                          ylim=c(min(Im(allVals)), max(Im(allVals))),
                          type="l", col="blue")
                  }
                }
              }
              
              par(mar=c(0,0,0,0))
              for (i in 1:p) {
                plot.new()
                text(0.5,0.5,substitute(paste(tau[1],"=",k),list(k=levels[i])), srt=90)
              }
              
              for (i in 1:p) {
                plot.new()
                text(0.5,0.5,substitute(paste(tau[2],"=",k),list(k=levels[i])))
              }
              if (length(xlab)>0) {
                plot.new()
                text(0.5, 0.5, xlab)
              }
              if (length(ylab)>0) {
                plot.new()
                text(0.5, 0.5, ylab, srt=90)
              }
              
            },  error = function(e) e,
            warning = function(w) w,
            finally = {
              par(def.par)  #- reset to default
            })
          }
)

  
  
