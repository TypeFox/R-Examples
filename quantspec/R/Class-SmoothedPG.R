#' @include generics.R
#' @include Class-Weight.R
NULL

################################################################################
#' Class for a smoothed quantile periodogram.
#'
#' \code{SmoothedPG} is an S4 class that implements the necessary
#' calculations to determine a smoothed version of one of the quantile
#' periodograms defined in Dette et. al (2015), Kley et. al (2016) and
#' Barunik&Kley (2015).
#'
#' For a \code{\link{QuantilePG}} \eqn{Q^{j_1, j_2}_n(\omega, x_1, x_2)}{Qn(w,x1,x2)} and
#' a \code{\link{Weight}} \eqn{W_n(\cdot)}{Wn(.)} the smoothed version
#' \deqn{\frac{2\pi}{n} \sum_{s=1}^{n-1} W_n(\omega-2\pi s / n) Q^{j_1, j_2}_n(2\pi s / n, x_1, x_2)}
#' is determined.
#'
#' The convolution required to determine the smoothed periodogram is implemented
#' using \code{\link[stats]{convolve}}.
#'
#' @name SmoothedPG-class
#' @aliases SmoothedPG
#' @exportClass SmoothedPG
#'
#' @keywords S4-classes
#'
#' @slot env An environment to allow for slots which need to be
#'       accessable in a call-by-reference manner:
#'       \describe{
#'         \item{\code{sdNaive}}{An array used for storage of the naively
#'             estimated standard deviations of the smoothed periodogram.}
#'         \item{\code{sdNaive.freq}}{a vector indicating for which frequencies
#' 							\code{sdNaive} has been computed so far.}
#'         \item{\code{sdNaive.done}}{a flag indicating whether \code{sdNaive}
#'             has been set yet.}
#'         \item{\code{sdBoot}}{An array used for storage of the standard
#'             deviations of the smoothed periodogram, estimated via
#'             bootstrap.}
#'         \item{\code{sdBoot.done}}{a flag indicating whether
#'             \code{sdBoot.naive} has been set yet.}
#'       }
#' @slot qPG the \code{\link{QuantilePG}} to be smoothed
#' @slot weight the \code{\link{Weight}} to be used for smoothing
#'
################################################################################
setClass(
    Class = "SmoothedPG",
    representation=representation(
        env = "environment",
        qPG = "QuantilePG",
        weight = "Weight"
    ),
    contains = "QSpecQuantity"
)

#' @importFrom stats convolve
setMethod(
    f = "initialize",
    signature = "SmoothedPG",
    definition = function(.Object, qPG, weight, frequencies, levels) {
      
      .Object@env <- new.env(parent=emptyenv())
      # .Object@env$sdNaive.done <- FALSE
      # .Object@env$sdCohNaive.done <- FALSE
      .Object@env$sdNaive.freq <- c()
      .Object@env$sdCohNaive.freq <- c()
      
      #.Object@env$sdCohNaive.freq.done <- c()
      
      .Object@env$sdBoot.done <- FALSE
      
      .Object@qPG <- qPG
      .Object@weight <- weight
      
      
      if (!(is.vector(frequencies)  && is.numeric(frequencies))) {
        stop("'frequencies' needs to be specified as a vector of real numbers")
      }
      
      N <- lenTS(.Object@qPG@freqRep@Y)
      
      .Object@levels <- levels
      K1 <- length(.Object@levels[[1]])
      K2 <- length(.Object@levels[[2]])
      
      # Transform all frequencies to Fourier frequencies
      # and remove redundancies
      if (class(weight) == "KernelWeight") {
        freq <- frequenciesValidator(frequencies, N, steps=1:6)
      } else if (class(weight) == "SpecDistrWeight") {
        freq <- frequenciesValidator(frequencies, N, steps=c(1:3,5:6))
      } else {
        stop("Cannot handle this type of Weight object.")
      }
      
      .Object@frequencies <- freq
      
      J <- length(freq)
      
#      N <- lenTS(Y)
#      K1 <- length(objLevels1)
#      K2 <- length(objLevels2)
#      D1 <- length(d1)
#      D2 <- length(d2)
      
      B <- .Object@qPG@freqRep@B
      P <- dim(.Object@qPG@freqRep@Y)[2]
      
      .Object@env$sdNaive <- array(0, dim=c(floor(N/2)+1, P, K1, P, K2) )
      .Object@env$sdCohNaive <- array(0, dim=c(floor(N/2)+1, P, K1, P, K2) ) # dim=c(floor(N/2)+1, D1, K1, D2, K2))
      .Object@env$sdCohSqNaive <- array(0, dim=c(floor(N/2)+1, P, K1, P, K2) )
      .Object@env$sigSq <- 0
      .Object@env$sdBoot <- array()
      
      .Object@values <- array(0, dim=c(J,P,K1,P,K2,B+1))
      
      if (class(weight) == "KernelWeight") {
        
        WW <- getValues(.Object@weight, N=N)
        Wnj <- weight@env$Wnj
        
        II <- getValues(.Object@qPG, levels.1 = levels[[1]], levels.2 = levels[[2]])
        II <- array(II, dim = c(N,P,K1,P,K2,B+1))
        
        if (max(freq) > 0) {
          
          # f performs convolution as defined in Brillinger (1975)
          f <- function(v) {
            A <- convolve(v[c(2:N,1)], rev(WW), type="c")[2:N]
            B <- WW[2:N] * v[1]
            return((A - B) / Wnj[1:(N-1)])
          }
          
          res <- (2*pi/N) * apply(II, c(2,3,4,5,6), f)
          
          posFreq <- freq[which(freq != 0)]
          pp <- round(N/(2*pi)*(posFreq %% (2*pi)))
          
          .Object@values[which(freq != 0),,,,,] <- res[pp,,,,,]
        }
        
        if (min(freq) == 0) {
          res <- apply(II, c(2,3,4,5,6), function(v) {sum(v[2:N]*rev(WW[2:N])) / Wnj[N]})
          .Object@values[which(freq == 0),,,,,] <- (2*pi/N) * res
        }
        rm(II)
      } else if (class(weight) == "SpecDistrWeight") {
        if (max(freq) > 0) {
          
          II <- getValues(.Object@qPG, frequencies = 2*pi*(1:(N-1))/N,
              levels.1 = levels[[1]], levels.2 = levels[[2]])
          #II <- array(II, dim = c(N,P,K1,P,K2,B+1))
          
          if (P == 1) {
            res <- (2*pi/N) * apply(II, c(2,3,4), cumsum)
            array(res, dim = c(N,P,K1,P,K2,B+1))
          } else {
            res <- (2*pi/N) * apply(II, c(2,3,4,5,6), cumsum)
          }
          
          posFreq <- freq[which(freq != 0)]
          pp <- round(N/(2*pi)*(posFreq %% (2*pi)))
          
          if (P == 1) {
            .Object@values[which(freq != 0),,,,,] <- res[pp,,,]
          } else {
            .Object@values[which(freq != 0),,,,,] <- res[pp,,,,,]
          }
        }
      }
      
      return(.Object)
    }
)

################################################################################
#' Get values from a smoothed quantile periodogram.
#'
#' The returned array of \code{values} is of dimension \code{[J,K1,K2,B+1]},
#' where \code{J=length(frequencies)}, \code{K1=length(levels.1)},
#' \code{K2=length(levels.2))}, and \code{B} denotes the
#' value stored in slot \code{B} of \code{freqRep} [that is the number of
#' boostrap repetitions performed on initialization].
#' At position \code{(j,k1,k2,b)}
#' the returned value is the one corresponding to \code{frequencies[j]},
#' \code{levels.1[k1]} and \code{levels.2[k2]} that are closest to the
#' \code{frequencies}, \code{levels.1} and \code{levels.2}
#' available in \code{object}; \code{\link{closest.pos}} is used to determine
#' what closest to means. \code{b==1} corresponds to the estimate without
#' bootstrapping; \code{b>1} corresponds to the \code{b-1}st bootstrap estimate.
#' 
#' If not only one, but multiple time series are under study, the dimension of
#' the returned vector is of dimension \code{[J,P,K1,P,K2,B+1]}, where \code{P}
#' denotes the dimension of the time series. 
#'
#' @name getValues-SmoothedPG
#' @aliases getValues,SmoothedPG-method
#'
#' @keywords Access-functions
#'
#' @param object \code{SmoothedPG} of which to get the values
#' @param frequencies a vector of frequencies for which to get the values
#' @param levels.1 the first vector of levels for which to get the values
#' @param levels.2 the second vector of levels for which to get the values
#' @param d1 optional parameter that determine for which j1 to return the
#' 					 data; may be a vector of elements 1, ..., D
#' @param d2 same as d1, but for j2
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
        "SmoothedPG"),
    definition = function(object,
        frequencies=2*pi*(0:(lenTS(object@qPG@freqRep@Y)-1))/lenTS(object@qPG@freqRep@Y),
        levels.1=getLevels(object,1),
        levels.2=getLevels(object,2),
        d1 = 1:(dim(object@values)[2]),
        d2 = 1:(dim(object@values)[4])) {
      
      # workaround: default values don't seem to work for generic functions?
      if (!hasArg(frequencies)) {
        frequencies <- 2*pi*(0:(lenTS(object@qPG@freqRep@Y)-1))/lenTS(object@qPG@freqRep@Y)
      }
      if (!hasArg(levels.1)) {
        levels.1 <- object@levels[[1]]
      }
      if (!hasArg(levels.2)) {
        levels.2 <- object@levels[[2]]
      }
      if (!hasArg(d1)) {
        d1 <- 1:(dim(object@values)[2])
      }
      if (!hasArg(d2)) {
        d2 <- 1:(dim(object@values)[4])
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
      #D <- dim(object@values)[2]
      D1 <- length(d1)
      D2 <- length(d2)
      K1 <- length(levels.1)
      K2 <- length(levels.2)
      res <- array(dim=c(J, D1, K1, D2, K2, object@qPG@freqRep@B+1))
      
      
      if (class(object@weight) == "KernelWeight") {
        
        # Select rows
        r1.pos <- closest.pos(oF, f[f <= pi])
        r2.pos <- closest.pos(-1*(2*pi-oF),-1*f[f > pi])
        
        if (length(r1.pos) > 0) {
          res[which(f <= pi),,,,,] <- object@values[r1.pos,d1,c.1.pos,d2,c.2.pos,]
        }
        if (length(r2.pos) > 0) {
          res[which(f > pi),,,,,] <- Conj(object@values[r2.pos,d1,c.1.pos,d2,c.2.pos,])
        }
        
      } else if (class(object@weight) == "SpecDistrWeight") {
        # Select rows
        r.pos <- closest.pos(oF, f)
        res[,,,,,] <- object@values[r.pos,,c.1.pos,,c.2.pos,]
      }
      
      if (D1 == 1 && D2 == 1) {
        final.dim.res <- c(J, K1, K2, object@qPG@freqRep@B+1)
      } else {
        final.dim.res <- c(J, D1, K1, D2, K2, object@qPG@freqRep@B+1)
      }
      
      res <- array(res, dim=final.dim.res)
      
      return(res)
    }
)

################################################################################
#' Compute quantile coherency from a smoothed quantile periodogram.
#'
#' Returns quantile coherency defined as
#' \deqn{\frac{G^{j_1, j_2}(\omega; \tau_1, \tau_2)}{(G^{j_1, j_1}(\omega; \tau_1, \tau_1) G^{j_2, j_2}(\omega; \tau_2, \tau_2))^{1/2}}}
#' where \eqn{G^{j_1, j_2}(\omega; \tau_1, \tau_2)} is the smoothed quantile
#' periodogram.
#' 
#' For the mechanism of selecting frequencies, dimensions and/or levels see,
#' for example, \code{\link{getValues-SmoothedPG}}.
#'
#' @name getCoherency-SmoothedPG
#' @aliases getCoherency,SmoothedPG-method
#'
#' @keywords Access-functions
#'
#' @param object \code{SmoothedPG} of which to get the values
#' @param frequencies a vector of frequencies for which to get the values
#' @param levels.1 the first vector of levels for which to get the values
#' @param levels.2 the second vector of levels for which to get the values
#' @param d1 optional parameter that determine for which j1 to return the
#' 					 data; may be a vector of elements 1, ..., D
#' @param d2 same as d1, but for j2
#'
#' @return Returns data from the array \code{values} that's a slot of
#'          \code{object}.
#'
#' @seealso
#' An example on how to use this function is analogously to the example given in
#' \code{\link{getValues-QuantilePG}}.
################################################################################
setMethod(f = "getCoherency",
    signature = signature(
        "SmoothedPG"),
    definition = function(object,
        frequencies=2*pi*(0:(lenTS(object@qPG@freqRep@Y)-1))/lenTS(object@qPG@freqRep@Y),
        levels.1=getLevels(object,1),
        levels.2=getLevels(object,2),
        d1 = 1:(dim(object@values)[2]),
        d2 = 1:(dim(object@values)[4])) {
      
      # workaround: default values don't seem to work for generic functions?
      if (!hasArg(frequencies)) {
        frequencies <- 2*pi*(0:(lenTS(object@qPG@freqRep@Y)-1))/lenTS(object@qPG@freqRep@Y)
      }
      if (!hasArg(levels.1)) {
        levels.1 <- object@levels[[1]]
      }
      if (!hasArg(levels.2)) {
        levels.2 <- object@levels[[2]]
      }
      if (!hasArg(d1)) {
        d1 <- 1:(dim(object@values)[2])
      }
      if (!hasArg(d2)) {
        d2 <- 1:(dim(object@values)[4])
      }
      # end: workaround
      
      J <- length(frequencies)
      #D <- dim(object@values)[2]
      D1 <- length(d1)
      D2 <- length(d2)
      K1 <- length(levels.1)
      K2 <- length(levels.2)
      res <- array(dim=c(J, D1, K1, D2, K2, object@qPG@freqRep@B+1))
      
      
      if (class(object@weight) != "KernelWeight") {
        stop("Coherency can only be determined if weight is of type KernelWeight.")
      }
      d <- union(d1,d2)
      V <- array(getValues(object, d1 = d, d2 = d, frequencies = frequencies), dim=c(J, D1, K1, D2, K2, object@qPG@freqRep@B+1))
      
      d1.pos <- closest.pos(d, d1)
      d2.pos <- closest.pos(d, d2)
      
      res <- .computeCoherency(V, d1.pos, d2.pos)
      
      if (D1 == 1 && D2 == 1) {
        final.dim.res <- c(J, K1, K2, object@qPG@freqRep@B+1)
      } else {
        final.dim.res <- c(J, D1, K1, D2, K2, object@qPG@freqRep@B+1)
      }

      res[is.nan(res)] <- NA
      res <- array(res, dim=final.dim.res)
      
      return(res)
    }
)

################################################################################
#' Get estimates for the standard deviation of the coherency computed from
#' smoothed quantile periodogram.
#'
#' Determines and returns an array of dimension \code{[J,K1,K2]},
#' where \code{J=length(frequencies)}, \code{K1=length(levels.1)}, and
#' \code{K2=length(levels.2))}. Whether
#' available or not, boostrap repetitions are ignored by this procedure.
#' At position \code{(j,k1,k2)}
#' the returned value is the standard deviation estimated corresponding to
#' \code{frequencies[j]}, \code{levels.1[k1]} and \code{levels.2[k2]} that are
#' closest to the
#' \code{frequencies}, \code{levels.1} and \code{levels.2}
#' available in \code{object}; \code{\link{closest.pos}} is used to determine
#' what closest to means.
#' 
#' If not only one, but multiple time series are under study, the dimension of
#' the returned vector is of dimension \code{[J,P,K1,P,K2]}, where \code{P}
#' denotes the dimension of the time series.
#'
#' Requires that the \code{\link{SmoothedPG}} is available at all Fourier
#' frequencies from \eqn{(0,\pi]}{(0,pi]}. If this is not the case the missing
#' values are imputed by taking one that is available and has a frequency
#' that is closest to the missing Fourier frequency; \code{closest.pos} is used
#' to determine which one this is.
#'
#' A precise definition on how the standard deviations of the smoothed quantile
#' periodogram are estimated is given in Barunik and Kley (2015). The estimate
#' returned is denoted by
#' \eqn{\sigma(\tau_1, \tau_2; \omega)}{sigma(tau1, tau2; omega)} on p. 26 of
#' the arXiv preprint.
#'
#' Note the ``standard deviation'' estimated here is not the square root of the
#' complex-valued variance. It's real part is the square root of the variance
#' of the real part of the estimator and the imaginary part is the square root
#' of the imaginary part of the variance of the estimator.
#'
#' @name getCoherencySdNaive-SmoothedPG
#' @aliases getCoherencySdNaive,SmoothedPG-method
#'
#' @keywords Access-functions
#'
#' @param object \code{\link{SmoothedPG}} of which to get the estimates for the
#'                standard deviation.
#' @param frequencies a vector of frequencies for which to get the result
#' @param levels.1 the first vector of levels for which to get the result
#' @param levels.2 the second vector of levels for which to get the result
#' @param d1 optional parameter that determine for which j1 to return the
#' 					 data; may be a vector of elements 1, ..., D
#' @param d2 same as d1, but for j2
#' @param type can be "1", where cov(Z, Conj(Z)) is subtracted, or "2", where
#' 					   it's not
#' @param impl choose "R" or "C" for one of the two implementations available
#'
#' @return Returns the estimate described above.
#'
#' @references
#' Kley, T., Volgushev, S., Dette, H. & Hallin, M. (2016).
#' Quantile Spectral Processes: Asymptotic Analysis and Inference.
#' \emph{Bernoulli}, \bold{22}(3), 1770--1807.
#' [cf. \url{http://arxiv.org/abs/1401.8104}]
#' 
#' Barunik, J. & Kley, T. (2015).
#' Quantile Cross-Spectral Measures of Dependence between Economic Variables.
#' [preprint available from the authors]
################################################################################
setMethod(f = "getCoherencySdNaive",
    signature = signature(
        object = "SmoothedPG"),
    definition = function(object,
        frequencies=2*pi*(0:(lenTS(object@qPG@freqRep@Y)-1))/lenTS(object@qPG@freqRep@Y),
        levels.1=getLevels(object,1),
        levels.2=getLevels(object,2),
        d1 = 1:(dim(object@values)[2]),
        d2 = 1:(dim(object@values)[4]),
        type = c("1", "2"),
        impl=c("R","C")) {
      
      if (class(getWeight(object)) != "KernelWeight") {
        stop("getSdNaive currently only available for 'KernelWeight'.")
      }
      
      Y <- object@qPG@freqRep@Y
      objLevels1 <- object@levels[[1]]
      objLevels2 <- object@levels[[2]]
      
      # workaround: default values don't seem to work for generic functions?
      if (!hasArg(frequencies)) {
        frequencies <- 2*pi*(0:(lenTS(Y)-1))/lenTS(Y)
      }
      if (!hasArg(levels.1)) {
        levels.1 <- objLevels1
      }
      if (!hasArg(levels.2)) {
        levels.2 <- objLevels2
      }
      if (!hasArg(d1)) {
        d1 <- 1:(dim(object@values)[2])
      }
      if (!hasArg(d2)) {
        d2 <- 1:(dim(object@values)[4])
      }
      if (!hasArg(type)) {
        type <- "1"
      }
      if (!hasArg(impl)) {
        impl <- "R"
      }
      # end: workaround
      
      N <- lenTS(Y)
      K1 <- length(objLevels1)
      K2 <- length(objLevels2)
      D1 <- length(d1)
      D2 <- length(d2)
      
      ## First determine which frequencies still have to be computed:
      #  i. e. for which j is 2*pi*j / N is to be computed.
      
      J.requested <- round(frequenciesValidator(frequencies, N = N) * N / (2*pi))  ## TODO: check if rounding required??
      J.toCompute <- setdiff(J.requested, object@env$sdCohNaive.freq)
      
      J <- length(J.toCompute)
      
      if ( J > 0 ) {
        
        # TODO: Make this work with type = 1 and type = 2
        #if (object@env$sdCohNaive.done == FALSE) {
        
        weight <- object@weight
        
        # Define a list which covariances to estimate
        # List shall be a matrix with rows
        # (ja ta jb tb jc tc jd td)
        #
        # by default all combinations shall be computed
        
        if (impl == "R") {
          
          lC <- matrix(ncol = 8, nrow = 7 * K1*K2*D1*D2 )
          
          i <- 1
          for (k1 in 1:K1) {
            for (i1 in 1:D1) {
              for (k2 in 1:K2) {
                for (i2 in 1:D2) {
                  jt_1 <- c(i1,k1)
                  jt_2 <- c(i2,k2)
                  lC[i,] <- c(jt_1, jt_1, jt_1, jt_1)
                  i <- i + 1
                  lC[i,] <- c(jt_1, jt_1, jt_1, jt_2)
                  i <- i + 1
                  lC[i,] <- c(jt_1, jt_1, jt_2, jt_2)
                  i <- i + 1
                  lC[i,] <- c(jt_1, jt_2, jt_1, jt_2)
                  i <- i + 1
                  lC[i,] <- c(jt_1, jt_2, jt_2, jt_1)
                  i <- i + 1
                  lC[i,] <- c(jt_1, jt_2, jt_2, jt_2)
                  i <- i + 1
                  lC[i,] <- c(jt_2, jt_2, jt_2, jt_2)
                  i <- i + 1
                }
              }
            }
          }
          
          lC <- matrix(lC[1:(i-1),], ncol = 8, nrow = i-1)
          lC <- unique(lC)
          
          #####
          ## Variant 1: more or less vectorized... 
          #####
          WW <- getValues(weight, N = N)[c(2:N,1)] # WW[j] corresponds to W_n(2 pi j/n)
          WW3 <- rep(WW,4) 
          
          # TODO: check whether it works for all d1, d2!!
          V <- array(getValues(object, frequencies = 2*pi*(1:(N-1))/N, d1=d1, d2=d2), dim=c(N-1,D1,K1,D2,K2))     
          res_coh <- array(0,dim=c(J,D1,K1,D2,K2))
          res_cohSq <- array(0,dim=c(J,D1,K1,D2,K2))
          auxRes <- array(0,dim=c(nrow(lC), J))
          
          M1 <- matrix(0, ncol=N-1, nrow=N)
          M2 <- matrix(0, ncol=N-1, nrow=N)
          #M3 <- matrix(0, ncol=N-1, nrow=N)
          #M4 <- matrix(0, ncol=N-1, nrow=N)
          for (j in 0:(N-1)) { # Test 1:N
            M1[j+1,] <- WW3[2*N+j-(1:(N-1))]*WW3[2*N+j-(1:(N-1))]
            M2[j+1,] <- WW3[2*N+j-(1:(N-1))]*WW3[2*N+j+(1:(N-1))]
            #M3[j+1,] <- WW3[2*N+j-(1:(N-1))]*WW3[2*N-j-(1:(N-1))]
            #M4[j+1,] <- WW3[2*N+j-(1:(N-1))]*WW3[2*N-j+(1:(N-1))]
          }
          
          M1 <- M1[J.toCompute+1,]
          M2 <- M2[J.toCompute+1,]
          
          for (r in 1:nrow(lC)) {
            jt <- lC[r,]
            V1 <- matrix(V[,jt[1],jt[2],jt[5],jt[6]] * Conj(V[,jt[3],jt[4],jt[7],jt[8]]), ncol=1)
            V2 <- matrix(V[,jt[1],jt[2],jt[7],jt[8]] * Conj(V[,jt[3],jt[4],jt[5],jt[6]]), ncol=1)
            
            auxRes[r,] <- rowSums(M1 %*% V1) + rowSums(M2 %*% V2)
          }
          
          V <- array(getValues(object, frequencies = 2*pi*(J.toCompute)/N, d1=d1, d2=d2), dim=c(J,D1,K1,D2,K2))
          Coh <- array(getCoherency(object, frequencies = 2*pi*(J.toCompute)/N, d1=d1, d2=d2), dim=c(J,D1,K1,D2,K2))
          
          
          for (i1 in 1:D1) {
            for (k1 in 1:K1) {
              for (i2 in 1:D2) {
                for (k2 in 1:K2) {
                  #r <- which(lC[,1] == i1 & lC[,2] == k1 & lC[,3] == i2 & lC[,4] == k2 & lC[,5] == i1 & lC[,6] == k1 & lC[,7] == i2 & lC[,8] == k2)
                  if (i1 == i2 && k1 == k2) {
                    S_coh <- 0
                    S_cohSq <- 0
                  } else {
                    H1111 <- auxRes[which(lC[,1] == i1 & lC[,2] == k1 & lC[,3] == i1 & lC[,4] == k1 & lC[,5] == i1 & lC[,6] == k1 & lC[,7] == i1 & lC[,8] == k1),]
                    H1112 <- auxRes[which(lC[,1] == i1 & lC[,2] == k1 & lC[,3] == i1 & lC[,4] == k1 & lC[,5] == i1 & lC[,6] == k1 & lC[,7] == i2 & lC[,8] == k2),]
                    H1122 <- auxRes[which(lC[,1] == i1 & lC[,2] == k1 & lC[,3] == i1 & lC[,4] == k1 & lC[,5] == i2 & lC[,6] == k2 & lC[,7] == i2 & lC[,8] == k2),]
                    H1212 <- auxRes[which(lC[,1] == i1 & lC[,2] == k1 & lC[,3] == i2 & lC[,4] == k2 & lC[,5] == i1 & lC[,6] == k1 & lC[,7] == i2 & lC[,8] == k2),]
                    H1221 <- auxRes[which(lC[,1] == i1 & lC[,2] == k1 & lC[,3] == i2 & lC[,4] == k2 & lC[,5] == i2 & lC[,6] == k2 & lC[,7] == i1 & lC[,8] == k1),]
                    H1222 <- auxRes[which(lC[,1] == i1 & lC[,2] == k1 & lC[,3] == i2 & lC[,4] == k2 & lC[,5] == i2 & lC[,6] == k2 & lC[,7] == i2 & lC[,8] == k2),]
                    H2222 <- auxRes[which(lC[,1] == i2 & lC[,2] == k2 & lC[,3] == i2 & lC[,4] == k2 & lC[,5] == i2 & lC[,6] == k2 & lC[,7] == i2 & lC[,8] == k2),]
                    
                    f12 <- V[,i1,k1,i2,k2]
                    f11 <- V[,i1,k1,i1,k1]
                    f22 <- V[,i2,k2,i2,k2]
                    
                    A <- H1212 - Re(f12*H1112/f11) - Re(f12*Conj(H1222)/f22)
                    A <- A + (abs(f12)^2/4) * (H1111/f11^2 + 2*Re(H1122/(f11*f22)) + H2222/f22^2)
                    A <- A/(f11*f22)
                    
                    B <- H1221 - f12*H1222/f22 - f12*Conj(H1112)/f11
                    B <- B + (f12^2/4) * (H1111/f11^2 + 2*Re(H1122/(f11*f22)) + H2222/f22^2)
                    B <- B/(f11*f22)
                    
                    #if (type == "1") {
                    S_coh <- (1/2) * complex(real = Re(A + B), imaginary = Re(A - B))
                    #} else {
                    # CORRECT??
                    # Recall, A == L1212 and B == L1221 ??
                    #R12 <- f12 / sqrt(abs(f11) * abs(f22)) # note that f11, f22 can be negative if weights are negative
                    R12 <- Coh[,i1,k1,i2,k2]
                    S_cohSq <- 2 * abs(R12)^2 * A
                    S_cohSq <- S_cohSq + 2 * (Re(R12)^2 - Im(R12)^2) * Re(B)
                    S_cohSq <- S_cohSq + 4 * Re(R12) * Im(R12) * Im(B)
                    #}
                    
                  }
                  ## TODO: Comment on the next line!!
                  ## Is this because S is always a linear combination of the
                  ## Cov(Lab, Lcd) terms??
                  
                  res_coh[,i1,k1,i2,k2] <- (2*pi/N)^2 * S_coh / ((weight@env$Wnj[c(N,1:(N-1))])[J.toCompute+1])^2
                  res_coh[,i2,k2,i1,k1] <- res_coh[,i1,k1,i2,k2]
                  
                  res_cohSq[,i1,k1,i2,k2] <- (2*pi/N)^2 * S_cohSq / ((weight@env$Wnj[c(N,1:(N-1))])[J.toCompute+1])^2
                  res_cohSq[,i2,k2,i1,k1] <- res_cohSq[,i1,k1,i2,k2]
                }
              }
            }
          }
          
          sqrt.cw <- function(z) {
            return(complex(real=sqrt(max(Re(z),1e-9)), imaginary=sqrt(max(Im(z),1e-9))))
          }
          #if (type == "1") {
          res_coh <- array(apply(res_coh,c(1,2,3,4,5),sqrt.cw), dim = c(J, D1, K1, D2, K2))  
          #} else {
          res_cohSq <- array(apply(res_cohSq,c(1,2,3,4,5),sqrt), dim = c(J, D1, K1, D2, K2))
          #}
          
          
          #####
          ## END Variant 1: more or less vectorized... 
          #####
          
        }
        
       
        object@env$sdCohNaive.freq <- union(object@env$sdCohNaive.freq, J.toCompute)
        object@env$sdCohNaive[J.toCompute+1,,,,] <- res_coh
        object@env$sdCohSqNaive[J.toCompute+1,,,,] <- res_cohSq
        
      } # End of if (J > 0)
#      } # End of 'if (object@env$sdNaive.done == FALSE) {'
      
      if (type == "1") {
        resObj <- object@env$sdCohNaive
      } else {
        resObj <- object@env$sdCohSqNaive
      }
      
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
      res <- array(dim=c(J, D1, K1, D2, K2))
      
      if (length(r1.pos) > 0) {
        res[which(f <= pi),,,,] <- resObj[r1.pos, , c.1.pos, , c.2.pos]
      }
      if (length(r2.pos) > 0) {
        res[which(f > pi),,,,] <- Conj(resObj[r2.pos, , c.1.pos, , c.2.pos])
      }
      
      return(res)
    }
)

################################################################################
#' Get estimates for the standard deviation of the smoothed quantile
#' periodogram.
#'
#' Determines and returns an array of dimension \code{[J,K1,K2]},
#' where \code{J=length(frequencies)}, \code{K1=length(levels.1)}, and
#' \code{K2=length(levels.2))}. Whether
#' available or not, boostrap repetitions are ignored by this procedure.
#' At position \code{(j,k1,k2)}
#' the returned value is the standard deviation estimated corresponding to
#' \code{frequencies[j]}, \code{levels.1[k1]} and \code{levels.2[k2]} that are
#' closest to the
#' \code{frequencies}, \code{levels.1} and \code{levels.2}
#' available in \code{object}; \code{\link{closest.pos}} is used to determine
#' what closest to means.
#' 
#' If not only one, but multiple time series are under study, the dimension of
#' the returned vector is of dimension \code{[J,P,K1,P,K2,B+1]}, where \code{P}
#' denotes the dimension of the time series. 
#'
#' Requires that the \code{\link{SmoothedPG}} is available at all Fourier
#' frequencies from \eqn{(0,\pi]}{(0,pi]}. If this is not the case the missing
#' values are imputed by taking one that is available and has a frequency
#' that is closest to the missing Fourier frequency; \code{closest.pos} is used
#' to determine which one this is.
#'
#' A precise definition on how the standard deviations of the smoothed quantile
#' periodogram are estimated is given in Barunik&Kley (2015).
#'
#' Note the ``standard deviation'' estimated here is not the square root of the
#' complex-valued variance. It's real part is the square root of the variance
#' of the real part of the estimator and the imaginary part is the square root
#' of the imaginary part of the variance of the estimator.
#'
#' @name getSdNaive-SmoothedPG
#' @aliases getSdNaive,SmoothedPG-method
#'
#' @keywords Access-functions
#'
#' @param object \code{\link{SmoothedPG}} of which to get the estimates for the
#'                standard deviation.
#' @param frequencies a vector of frequencies for which to get the result
#' @param levels.1 the first vector of levels for which to get the result
#' @param levels.2 the second vector of levels for which to get the result
#' @param d1 optional parameter that determine for which j1 to return the
#' 					 data; may be a vector of elements 1, ..., D
#' @param d2 same as d1, but for j2
#' @param impl choose "R" or "C" for one of the two implementations available
#'
#' @return Returns the estimate described above.
#'
#' @references
#' Dette, H., Hallin, M., Kley, T. & Volgushev, S. (2015).
#' Of Copulas, Quantiles, Ranks and Spectra: an \eqn{L_1}{L1}-approach to
#' spectral analysis. \emph{Bernoulli}, \bold{21}(2), 781--831.
#' [cf. \url{http://arxiv.org/abs/1111.7205}]
################################################################################
setMethod(f = "getSdNaive",
    signature = signature(
        object = "SmoothedPG"),
    definition = function(object,
        frequencies=2*pi*(0:(lenTS(object@qPG@freqRep@Y)-1))/lenTS(object@qPG@freqRep@Y),
        levels.1=getLevels(object,1),
        levels.2=getLevels(object,2),
        d1 = 1:(dim(object@values)[2]),
        d2 = 1:(dim(object@values)[4]),
        impl=c("R","C")) {
      
      if (class(getWeight(object)) != "KernelWeight") {
        stop("getSdNaive currently only available for 'KernelWeight'.")
      }
      
      Y <- object@qPG@freqRep@Y
      objLevels1 <- object@levels[[1]]
      objLevels2 <- object@levels[[2]]
      
      # workaround: default values don't seem to work for generic functions?
      if (!hasArg(frequencies)) {
        frequencies <- 2*pi*(0:(lenTS(Y)-1))/lenTS(Y)
      }
      if (!hasArg(levels.1)) {
        levels.1 <- objLevels1
      }
      if (!hasArg(levels.2)) {
        levels.2 <- objLevels2
      }
      if (!hasArg(d1)) {
        d1 <- 1:(dim(object@values)[2])
      }
      if (!hasArg(d2)) {
        d2 <- 1:(dim(object@values)[4])
      }
      if (!hasArg(impl)) {
        impl <- "R"
      }
      # end: workaround
      
      N <- lenTS(Y)
      
      #TODO: modify this code such that only the levels and dimensions
      # 		 that were not computed before are included in the computation.
      K1 <- length(objLevels1)
      K2 <- length(objLevels2)
      D1 <- ncol(Y)
      D2 <- ncol(Y)
      
      ## First determine which frequencies still have to be computed:
      #  i. e. for which j is 2*pi*k / N is to be computed.
      
      J.requested <- round(frequenciesValidator(frequencies, N = N) * N / (2*pi))  ## TODO: check if rounding required??
      J.toCompute <- setdiff(J.requested, object@env$sdNaive.freq)
      
      J <- length(J.toCompute)
      
      if ( J > 0 ) {
        #if (object@env$sdNaive.done == FALSE) {
        #if (!identical(object@env$sdNaive.freq.done, frequencies)) {
        
        weight <- object@weight
        
        # Define a list which covariances to estimate
        # List shall be a matrix with rows
        # (ja ta jb tb jc tc jd td)
        #
        # by default all combinations shall be computed
        
        if (impl == "R") {
          
          lC <- matrix(ncol = 8, nrow = K1*K2*D1*D2 + K1*D1*(K2*D2-1))
          
          i <- 1
          for (k1 in 1:K1) {
            for (i1 in 1:D1) {
              for (k2 in 1:K2) {
                for (i2 in 1:D2) {
                  lC[i,] <- rep(c(i1,k1,i2,k2),2)
                  i <- i + 1
                  if (i1 != i2 || k1 != k2) {
                    lC[i,] <- c(i1,k1,i2,k2,i2,k2,i1,k1)
                    i <- i + 1 
                  }
                }
              }
            }
          }
          
          #lC <- matrix(lC[1:(i-1),], ncol = 8, nrow = i-1)
          
          
          #####
          ## Variant 1: more or less vectorized... 
          #####
          WW <- getValues(weight, N = N)[c(2:N,1)] # WW[j] corresponds to W_n(2 pi j/n)
          WW3 <- rep(WW,4) 
          
          # TODO: fix to make it work for all d1, d2!!
          V <- array(getValues(object, frequencies = 2*pi*(1:(N-1))/N), dim=c(N-1,D1,K1,D2,K2))     
          res <- array(0,dim=c(J,D1,K1,D2,K2))
          auxRes <- array(0,dim=c(K1*K2*D1*D2 + K1*D1*(K2*D2-1), J))
          
          M1 <- matrix(0, ncol=N-1, nrow=N)
          M2 <- matrix(0, ncol=N-1, nrow=N)
          #M3 <- matrix(0, ncol=N-1, nrow=N)
          #M4 <- matrix(0, ncol=N-1, nrow=N)
          for (j in 0:(N-1)) { # Test 1:N
            M1[j+1,] <- WW3[2*N+j-(1:(N-1))]*WW3[2*N+j-(1:(N-1))]
            M2[j+1,] <- WW3[2*N+j-(1:(N-1))]*WW3[2*N+j+(1:(N-1))]
            #M3[j+1,] <- WW3[2*N+j-(1:(N-1))]*WW3[2*N-j-(1:(N-1))]
            #M4[j+1,] <- WW3[2*N+j-(1:(N-1))]*WW3[2*N-j+(1:(N-1))]
          }
          
          M1 <- M1[J.toCompute+1,]
          M2 <- M2[J.toCompute+1,]
          
          for (r in 1:nrow(lC)) {
            jt <- lC[r,]
            V1 <- matrix(V[,jt[1],jt[2],jt[5],jt[6]] * Conj(V[,jt[3],jt[4],jt[7],jt[8]]), ncol=1)
            V2 <- matrix(V[,jt[1],jt[2],jt[7],jt[8]] * Conj(V[,jt[3],jt[4],jt[5],jt[6]]), ncol=1)
            
            auxRes[r,] <- rowSums(M1 %*% V1) + rowSums(M2 %*% V2)
            #auxRes[r,2,] <- rowSums(M3 %*% V1) + rowSums(M4 %*% V2)
          }
          
          for (i1 in 1:D1) {
            for (k1 in 1:K1) {
              for (i2 in 1:D2) {
                for (k2 in 1:K2) {
                  r <- which(lC[,1] == i1 & lC[,2] == k1 & lC[,3] == i2 & lC[,4] == k2 & lC[,5] == i1 & lC[,6] == k1 & lC[,7] == i2 & lC[,8] == k2)
                  if (i1 == i2 && k1 == k2) {
                    S <- auxRes[r,]
                  } else {
                    r2 <- which(lC[,1] == i1 & lC[,2] == k1 & lC[,3] == i2 & lC[,4] == k2 & lC[,5] == i2 & lC[,6] == k2 & lC[,7] == i1 & lC[,8] == k1)
                    S <- (1/2) * complex(real = Re(auxRes[r,] + auxRes[r2,]), imaginary = Re(auxRes[r,] - auxRes[r2,]))
                  }
                  
                  #res[,i1,k1,i2,k2] <- (2*pi/N)^2 * S / (weight@env$Wnj[c(N,1:(N-1))])^2
                  res[,i1,k1,i2,k2] <- (2*pi/N)^2 * S / ((weight@env$Wnj[c(N,1:(N-1))])[J.toCompute+1])^2
                  res[,i2,k2,i1,k1] <- res[,i1,k1,i2,k2]
                }
              }
            }
          }
          
          sqrt.cw <- function(z) {
            return(complex(real=sqrt(max(Re(z),0)), imaginary=sqrt(max(Im(z),0))))
          }
          res <- array(apply(res,c(1,2,3,4,5),sqrt.cw), dim = c(J, D1, K1, D2, K2))
          
          #####
          ## END Variant 1: more or less vectorized... 
          #####
          
        }
        
#        if (impl == "C") {     
#          #####
#          ## Variant 2: using C++
#          #####
#          
#          WW <- getValues(weight, N=N) 
#          V <- array(getValues(object, frequencies = 2*pi*(0:(N-1))/N)[,,,,,1], dim=c(N,K1,K2))     
#          res <- .computeSdNaive(V, WW)
#     
#          #####
#          ## END Variant 2: using C++ (cppFunction)
#          #####
#          
#        }
        
        object@env$sdNaive.freq <- union(object@env$sdNaive.freq, J.toCompute)
        
        object@env$sdNaive[J.toCompute+1,,,,] <- res
        
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
      res <- array(dim=c(J, length(d1), K1, length(d2), K2))
      
      if (length(r1.pos) > 0) {
        res[which(f <= pi),,,,] <- resObj[r1.pos, d1, c.1.pos, d2, c.2.pos]
      }
      if (length(r2.pos) > 0) {
        res[which(f > pi),,,,] <- Conj(resObj[r2.pos, d1, c.1.pos, d2, c.2.pos])
      }    
      
      if (length(d1) == 1 && length(d2) == 1) {
        final.dim.res <- c(J, K1, K2)
      } else {
        final.dim.res <- c(J, length(d1), K1, length(d2), K2)
      }
      
      res <- array(res, dim=final.dim.res)
      
      return(res)
    }
)

################################################################################
#' Get bootstrap estimates for the standard deviation of the smoothed quantile
#' periodogram.
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
#' Requires that the \code{\link{SmoothedPG}} is available at all Fourier
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
#' @name getSdBoot-SmoothedPG
#' @aliases getSdBoot,SmoothedPG-method
#' 
#' @importFrom stats sd
#'
#' @keywords Access-functions
#'
#' @param object \code{\link{SmoothedPG}} of which to get the bootstrap estimates for the
#'                standard deviation.
#' @param frequencies a vector of frequencies for which to get the result
#' @param levels.1 the first vector of levels for which to get the result
#' @param levels.2 the second vector of levels for which to get the result
#'
#' @return Returns the estimate described above.
################################################################################

setMethod(f = "getSdBoot",
    signature = signature(
        object = "SmoothedPG"),
    definition = function(object,
        frequencies=2*pi*(0:(lenTS(object@qPG@freqRep@Y)-1))/lenTS(object@qPG@freqRep@Y),
        levels.1=getLevels(object,1),
        levels.2=getLevels(object,2)) {
      
      if (class(getWeight(object)) != "KernelWeight") {
        stop("getSdBoot currently only available for 'KernelWeight'.")
      }
      
      # workaround: default values don't seem to work for generic functions?
      if (!hasArg(frequencies)) {
        frequencies <- 2*pi*(0:(lenTS(object@qPG@freqRep@Y)-1))/lenTS(object@qPG@freqRep@Y)
      }
      if (!hasArg(levels.1)) {
        levels.1 <- object@levels[[1]]
      }
      if (!hasArg(levels.2)) {
        levels.2 <- object@levels[[2]]
      }
      # end: workaround
      
      if (dim(object@env$sdBoot) == 1) {
        
        complex.var <- function(x) {
          return(complex(real = sd(Re(x)), imaginary = sd(Im(x))))
        }
        
        #N <- length(object@qPG@freqRep@Y)
        #K1 <- length(object@levels[[1]])
        #K2 <- length(object@levels[[2]])
        B <- object@qPG@freqRep@B
        
        # TODO: fix...
        v <- getValues(object, frequencies = frequencies,
            levels.1 = levels.1, levels.2 = levels.2, d1=1, d2=1)[,,,2:(B+1), drop=F]
        object@env$sdBoot <- apply(v, c(1,2,3), complex.var)
      }
      return(object@env$sdBoot)
    }
)

################################################################################
#' Get pointwise confidence intervals for the quantile spectral density kernel,
#' quantile coherency or quantile coherence.
#'
#' Returns a list of two arrays \code{lowerCIs} and \code{upperCIs} that contain
#' the upper and lower limits for a level \code{1-alpha} confidence interval of
#' the quantity of interest. Each array is of dimension \code{[J,K1,K2]} if a
#' univariate time series is being analysed or of dimension \code{[J,D1,K1,D2,K2]},
#' where \code{J=length(frequencies)}, \code{D1=length(d1)}, \code{D2=length(d2)},
#' \code{K1=length(levels.1)}, and \code{K2=length(levels.2))}.
#' At position \code{(j,k1,k2)} or \code{(j,i1,k1,i2,k2)} the real (imaginary)
#' part of the returned values are the bounds of the confidence interval for the
#' the real (imaginary) part of the quantity under anlysis, which corresponds to
#' \code{frequencies[j]}, \code{d1[i1]}, \code{d2[i2]}, \code{levels.1[k1]} and
#' \code{levels.2[k2]} closest to the Fourier frequencies, \code{levels.1} and
#' \code{levels.2} available in \code{object}; \code{\link{closest.pos}} is used
#' to determine what closest to means.
#'
#' Currently, pointwise confidence bands for two different \code{quantity}
#' are implemented:
#' \itemize{
#'   \item \code{"spectral density"}: confidence intervals for the quantile spectral
#' 					 density as described in Kley et. al (2016) for the univariate case and
#' 					 in Barunik and Kley (2015) for the multivariate case.
#'   \item \code{"coherency"}: confidence intervals for the quantile coherency as
#' 					 described in Barunik and Kley (2015).
#' }
#' 
#' Currently, three different \code{type}s of confidence intervals are
#' available:
#' \itemize{
#'   \item \code{"naive.sd"}: confidence intervals based on the asymptotic
#'           normality of the smoothed quantile periodogram; standard deviations
#'           are estimated using \code{\link{getSdNaive}}.
#'   \item \code{"boot.sd"}: confidence intervals based on the asymptotic
#'           normality of the smoothed quantile periodogram; standard deviations
#'           are estimated using \code{\link{getSdBoot}}.
#'   \item \code{"boot.full"}: confidence intervals determined by estimating the
#'           quantiles of he distribution of the smoothed quantile periodogram,
#'           by the empirical quantiles of the sample of bootstrapped
#'           replications.
#' }
#'
#' @name getPointwiseCIs-SmoothedPG
#' @aliases getPointwiseCIs,SmoothedPG-method
#' 
#' @importFrom stats qnorm
#' @importFrom stats quantile
#'
#' @keywords Access-functions
#'
#' @param object \code{SmoothedPG} of which to get the confidence intervals
#' @param quantity a flag indicating for which the pointwise confidence bands
#' 								 will be determined. Can take one of the possible values
#' 								 discussed above. 
#' @param frequencies a vector of frequencies for which to get the result
#' @param levels.1 the first vector of levels for which to get the result
#' @param levels.2 the second vector of levels for which to get the result
#' @param d1 optional parameter that determine for which j1 to return the
#' 					 data; may be a vector of elements 1, ..., D
#' @param d2 same as d1, but for j2
#' @param alpha the level of the confidence interval; must be from \eqn{(0,1)}
#' @param type a flag indicating which type of confidence interval should be
#'         returned; can take one of the three values discussed above.
#'
#' @return Returns a named list of two arrays \code{lowerCIS} and \code{upperCIs}
#'          containing the lower and upper bounds for the confidence intervals.
#'
#' @examples
#' sPG <- smoothedPG(rnorm(2^10), levels.1=0.5)
#' CI.upper <- Re(getPointwiseCIs(sPG)$upperCIs[,1,1])
#' CI.lower <- Re(getPointwiseCIs(sPG)$lowerCIs[,1,1])
#' freq = 2*pi*(0:1023)/1024
#' plot(x = freq, y = rep(0.25/(2*pi),1024),
#'    ylim=c(min(CI.lower), max(CI.upper)),
#'    type="l", col="red") # true spectrum
#' lines(x = freq, y = CI.upper)
#' lines(x = freq, y = CI.lower)
################################################################################
setMethod(f = "getPointwiseCIs",
    signature = signature(
        object = "SmoothedPG"),
    definition = function(object,
        quantity = c("spectral density", "coherency", "coherence"),
        frequencies=2*pi*(0:(lenTS(object@qPG@freqRep@Y)-1))/lenTS(object@qPG@freqRep@Y),
        levels.1=getLevels(object,1),
        levels.2=getLevels(object,2),
        d1 = 1:(dim(object@values)[2]),
        d2 = 1:(dim(object@values)[4]),
        alpha=.1, type=c("naive.sd", "boot.sd", "boot.full")) {
      
      # workaround: default values don't seem to work for generic functions?
      if (!hasArg(frequencies)) {
        frequencies <- 2*pi*(0:(lenTS(object@qPG@freqRep@Y)-1))/lenTS(object@qPG@freqRep@Y)
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
      if (!hasArg(d1)) {
        d1 <- 1:(dim(object@values)[2])
      }
      if (!hasArg(d2)) {
        d2 <- 1:(dim(object@values)[4])
      }
      if (!hasArg(type)) {
        type <- "naive.sd"
      }
      # end: workaround
      
      type <- match.arg(type)[1]
      quantity <- match.arg(quantity)[1]
      switch(type,
          "naive.sd" = {
            switch(quantity,
                "spectral density" = {
                    sdEstim <- getSdNaive(object,
                        frequencies = frequencies,
                        levels.1 = levels.1,
                        levels.2 = levels.2,
                        d1 = d1, d2 = d2)
                  },
                "coherency" = {
                    sdEstim <- getCoherencySdNaive(object,
                        frequencies = frequencies,
                        levels.1 = levels.1,
                        levels.2 = levels.2,
                        d1 = d1, d2 = d2, type="1")
                },
                "coherence" = {
                    sdEstim <- getCoherencySdNaive(object,
                        frequencies = frequencies,
                        levels.1 = levels.1,
                        levels.2 = levels.2,
                        d1 = d1, d2 = d2, type="2")
                })
            },
          "boot.sd" = {
            if (quantity == "spectral density") {
              sdEstim <- getSdBoot(object,
                  frequencies = frequencies,
                  levels.1 = levels.1,
                  levels.2 = levels.2)
            } else {
              stop("boot.sd is so far only implemented for quantity = 'spectral density'.")
            }
        }
      )
      
      J <- length(frequencies)
      K1 <- length(levels.1)
      K2 <- length(levels.2)
      D1 <- length(d1)
      D2 <- length(d2)
      
      
      if (type == "naive.sd" || type == "boot.sd") {
        
        switch(quantity,
            "spectral density" = {
              v <- array(getValues(object,
                      frequencies = frequencies,
                      levels.1 = levels.1,
                      levels.2 = levels.2, d1=d1, d2=d2), dim = c(J, D1, K1, D2, K2))
              sdEstim <- array(sdEstim, dim = c(J, D1, K1, D2, K2))
              upperCIs <- array(v + sdEstim * qnorm(1-alpha/2), dim = c(J, D1, K1, D2, K2))
              lowerCIs <- array(v + sdEstim * qnorm(alpha/2), dim = c(J, D1, K1, D2, K2))
            },
            "coherency" = {
              v <- array(getCoherency(object,
                      frequencies = frequencies,
                      levels.1 = levels.1,
                      levels.2 = levels.2, d1=d1, d2=d2), dim = c(J, D1, K1, D2, K2))
              upperCIs <- array(v + sdEstim * qnorm(1-alpha/2), dim = c(J, D1, K1, D2, K2))
              lowerCIs <- array(v + sdEstim * qnorm(alpha/2), dim = c(J, D1, K1, D2, K2))
            },
            "coherence" = {
              v <- array(getCoherency(object,
                      frequencies = frequencies,
                      levels.1 = levels.1,
                      levels.2 = levels.2, d1=d1, d2=d2), dim = c(J, D1, K1, D2, K2))
              
              upperCIs <- array(abs(v)^2 + Re(sdEstim) * qnorm(1-alpha/2), dim = c(J, D1, K1, D2, K2))
              lowerCIs <- array(abs(v)^2 + Re(sdEstim) * qnorm(alpha/2), dim = c(J, D1, K1, D2, K2))
            })

      } else if (type == "boot.full") {
        
        if (quantity == "spectral density") {
          
          # TODO: Error Msg ausgeben falls B == 0
          B <- object@qPG@freqRep@B
          # TODO: fix...
    
    switch(quantity,
        "spectral density" = {
          v <- getValues(object,
              frequencies = frequencies,
              levels.1 = levels.1,
              levels.2 = levels.2, d1=1, d2=1)[,,,2:(B+1), drop=FALSE]

        },
        "coherency" = {
          v <- getCoherency(object,
              frequencies = frequencies,
              levels.1 = levels.1,
              levels.2 = levels.2, d1=1, d2=1)[,,,2:(B+1), drop=FALSE]
        },
        "coherence" = {
          v <- getCoherency(object,
              frequencies = frequencies,
              levels.1 = levels.1,
              levels.2 = levels.2, d1=1, d2=1)[,,,2:(B+1), drop=FALSE]
          v <- abs(v)^2
        })
    
          v <- getValues(object,
              frequencies = frequencies,
              levels.1 = levels.1,
              levels.2 = levels.2, d1=1, d2=1)[,,,2:(B+1), drop=FALSE]
          uQuantile <- function(x) {complex(real = quantile(Re(x),1-alpha/2),
                imaginary = quantile(Im(x),1-alpha/2))}
          lQuantile <- function(x) {complex(real = quantile(Re(x),alpha/2),
                imaginary = quantile(Im(x),alpha/2))}
          
          upperCIs <- apply(v, c(1,2,3), uQuantile)
          lowerCIs <- apply(v, c(1,2,3), lQuantile)
          
        }
      }
      
      if (D1 == 1 && D2 == 1) {
        final.dim.res <- c(J, K1, K2)
      } else {
        final.dim.res <- c(J, D1, K1, D2, K2)
      }
      
      lowerCIs <- array(lowerCIs, dim=final.dim.res)
      upperCIs <- array(upperCIs, dim=final.dim.res)
      
      res <- list(lowerCIs = lowerCIs, upperCIs = upperCIs)
      return(res)

    }
)

################################################################################
#' Get associated \code{\link{Weight}} from a \code{\link{SmoothedPG}}.
#'
#' @name getWeight-SmoothedPG
#' @aliases getWeight,SmoothedPG-method
#'
#' @keywords Access-association-functions
#'
#' @param object \code{SmoothedPG} from which to get the \code{Weight}.
#' @return Returns the \code{\link{Weight}} object associated.
################################################################################
setMethod(f = "getWeight",
    signature = "SmoothedPG",
    definition = function(object) {
      return(object@weight)
    }
)

################################################################################
#' Get associated \code{\link{QuantilePG}} from a \code{\link{SmoothedPG}}.
#'
#' @name getQuantilePG-SmoothedPG
#' @aliases getQuantilePG,SmoothedPG-method
#'
#' @keywords Access-association-functions
#'
#' @param object \code{SmoothedPG} from which to get the \code{\link{QuantilePG}}.
#' @return Returns the \code{\link{QuantilePG}} object associated.
################################################################################
setMethod(f = "getQuantilePG",
    signature = "SmoothedPG",
    definition = function(object) {
      return(object@qPG)
    }
)

################################################################################
#' Create an instance of the \code{SmoothedPG} class.
#'
#' A \code{SmoothedPG} object can be created from either
#' \itemize{
#'  \item a \code{numeric}, a \code{ts}, or a \code{zoo} object
#'   \item a \code{QuantilePG} object.
#' }
#' If a \code{QuantilePG} object is used for smoothing, only the \code{weight},
#' \code{frequencies} and \code{levels.1} and \code{levels.2} parameters are
#' used; all others are ignored. In this case the default values for the levels
#' are the levels of the \code{QuantilePG} used for smoothing. Any subset of the
#' levels available there can be chosen.
#'
#' The parameter \code{type.boot} can be set to choose a block bootstrapping
#' procedure. If \code{"none"} is chosen, a moving blocks bootstrap with
#' \code{l=length(Y)} and  \code{N=length(Y)} would be done. Note that in that
#' case one would also chose \code{B=0} which means that \code{getPositions}
#' would never be called. If \code{B>0} then each bootstrap replication would
#' be the undisturbed time series.
#'
#' @name SmoothedPG-constructor
#' @aliases smoothedPG
#' @export
#'
#' @keywords Constructors
#'
#' @param object a time series (\code{numeric}, \code{ts}, or \code{zoo} object)
#'                from which to determine the smoothed periodogram; alternatively
#'                a \code{\link{QuantilePG}} object can be supplied.
#' @param isRankBased If true the time series is first transformed to pseudo
#'                    data [cf. \code{\link{FreqRep}}].
#' @param levels.1 A vector of length \code{K1} containing the levels \code{x1}
#'                  at which the SmoothedPG is to be determined.
#' @param levels.2 A vector of length \code{K2} containing the levels \code{x2}.
#' @param frequencies A vector containing frequencies at which to determine the
#'                     smoothed periodogram.
#' @param type A flag to choose the type of the estimator. Can be either
#'              \code{"clipped"} or \code{"qr"}. In the first case
#'              \code{\link{ClippedFT}} is used as a frequency representation, in
#'              the second case \code{\link{QRegEstimator}} is used.
#' @param B number of bootstrap replications
#' @param l (expected) length of blocks
#' @param type.boot A flag to choose a method for the block bootstrap; currently
#'                   two options are implemented: \code{"none"} and \code{"mbb"}
#'                   which means to do a moving blocks  bootstrap with \code{B}
#'                   and \code{l} as specified.
#' @param method  method used for computing the quantile regression estimates.
#'                 The choice is passed to \code{qr}; see the
#'                 documentation of \code{quantreg} for details.
#' @param parallel a flag to allow performing parallel computations,
#'                   where possible.
#' @param weight Object of type \code{\link{Weight}} to be used for smoothing.
#'
#' @return Returns an instance of \code{SmoothedPG}.
#'
#' @examples
#' Y <- rnorm(64)
#' levels.1 <- c(0.25,0.5,0.75)
#' weight <- kernelWeight(W=W0)
#'
#' # Version 1a of the constructor -- for numerics:
#' sPG.ft <- smoothedPG(Y, levels.1 = levels.1, weight = weight, type="clipped")
#' sPG.qr <- smoothedPG(Y, levels.1 = levels.1, weight = weight, type="qr")
#'
#' # Version 1b of the constructor -- for ts objects:
#' sPG.ft <- smoothedPG(wheatprices, levels.1 = c(0.05,0.5,0.95), weight = weight)
#'
#' # Version 1c of the constructor -- for zoo objects:
#' sPG.ft <- smoothedPG(sp500, levels.1 = c(0.05,0.5,0.95), weight = weight)
#'
#' # Version 2 of the constructor:
#' qPG.ft <- quantilePG(Y, levels.1 = levels.1, type="clipped")
#' sPG.ft <- smoothedPG(qPG.ft, weight = weight)
#' qPG.qr <- quantilePG(Y, levels.1 = levels.1, type="qr")
#' sPG.qr <- smoothedPG(qPG.qr, weight = weight)
################################################################################
smoothedPG <- function(
    object,
    frequencies=2*pi/lenTS(object) * 0:(lenTS(object)-1),
    levels.1 = 0.5,
    levels.2=levels.1,
    isRankBased=TRUE,
    type=c("clipped","qr"),
    type.boot=c("none","mbb"),
    method = c("br", "fn", "pfn", "fnc", "lasso", "scad"),
    parallel=FALSE,
    B = 0,
    l = 1,
    weight = kernelWeight()) {
  
  if (class(object) == "numeric" | class(object) == "matrix") {
    versConstr <- 1
    Y <- object
  } else if (class(object) == "ts") {
    versConstr <- 1
    Y <- object[1:length(object)]
  } else if (class(object) == "zoo") {
    versConstr <- 1
    Y <- coredata(object)
  } else if (class(object) == "QuantilePG") {
    versConstr <- 2
    
    
    
    if (!hasArg(frequencies)) {
      Y <- object@freqRep@Y
      frequencies <- 2*pi/lenTS(Y) * 0:(lenTS(Y)-1)
    }
    
    if (!hasArg(levels.1)) {
      levels.1 <- getLevels(object,1)
    }
    if (!hasArg(levels.2)) {
      levels.2 <- getLevels(object,2)
    }
    
  } else {
    stop("object is neither 'numeric', 'matrix', 'ts', 'zoo', nor 'QuantilePG'.")
  }
  
  if (versConstr == 1) {
    
    levels.all <- union(levels.1, levels.2)
    K1 <- length(levels.1)
    K2 <- length(levels.2)
    
    J <- length(frequencies)
    
    qPG <- quantilePG(Y, frequencies = 2*pi/lenTS(Y) * 0:(lenTS(Y)-1),
        levels.1 = levels.1,
        levels.2 = levels.2,
        isRankBased = isRankBased,
        type=type, type.boot = type.boot, B=B, l=l)
  } else if (versConstr == 2) {
    qPG <- object
  }
  
  obj <- new(
      Class = "SmoothedPG",
      qPG = qPG,
      weight = weight,
      frequencies = frequencies,
      levels = list(levels.1, levels.2)
  )
  
  return(obj)
  
}


################################################################################
#' Plot the values of a \code{\link{SmoothedPG}}.
#'
#' Creates a \code{K} x \code{K} plot depicting a smoothed quantile periodogram.
#' Optionally, the quantile periodogram on which the smoothing was performed,
#' a simulated quantile spectral density, and pointwise confidence intervals can
#' be displayed.
#' In each of the subplots either the real part (on and below the diagonal;
#' i. e., \eqn{\tau_1 \leq \tau_2}{tau1 <= tau2}) or the imaginary parts
#' (above the diagonal; i. e., \eqn{\tau_1 > \tau_2}{tau1 > tau2}) of
#' \itemize{
#'   \item the smoothed quantile periodogram (blue line),
#'   \item the quanitle peridogram that was used for smoothing (gray line),
#'   \item a simulated quantile spectral density (red line),
#'   \item pointwise (asymptotic) confidence intervals (light gray area),
#' }
#' for the combination of levels \eqn{\tau_1}{tau1} and \eqn{\tau_2}{tau2}
#' denoted on the left and bottom margin of the plot are displayed.
#' 
#' Currently, only the plot for the first component is shown.
#'
#' @name plot-SmoothedPG
#' @aliases plot,SmoothedPG,ANY-method
#' @export
#'
#' @importFrom abind abind
#' @importFrom grDevices gray
#'
#' @param x  The \code{\link{SmoothedPG}} object to plot
#' @param plotPG a flag indicating weater the \code{QuantilePG} object
#'                associated with the \code{\link{SmoothedPG}} \code{x}
#'                is also to be plotted.
#' @param qsd  a \code{\link{QuantileSD}} object; will be plotted if not
#'              missing.
#' @param ptw.CIs the confidence level for the confidence intervals to be
#'                 displayed; must be a number from [0,1]; if null, then no
#'                 confidence intervals will be plotted.
#' @param type.CIs indicates the method to be used for determining the
#'                 confidence intervals; the methods available are those
#'                  provided by
#'                 \code{\link{getPointwiseCIs-SmoothedPG}}.
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
#' @param frequencies a set of frequencies for which the values are to be
#'                    plotted.
#' @param levels a set of levels for which the values are to be plotted.
#'
#' @return Returns the plot described in the Description section.
################################################################################

setMethod(f = "plot",
    signature = signature("SmoothedPG"),
    definition = function(x, plotPG=FALSE, qsd,
        ptw.CIs = 0.1, type.CIs = c("naive.sd", "boot.sd", "boot.full"),
        ratio = 3/2, widthlab = lcm(1), xlab = expression(omega/2*pi), ylab = NULL,
        type.scaling = c("individual", "real-imaginary", "all"),
        frequencies=x@frequencies,
        levels=intersect(x@levels[[1]], x@levels[[2]])) {
      
      def.par <- par(no.readonly = TRUE) # save default, for resetting...
      
      # workaround: default values don't seem to work for generic functions?
      if (!hasArg(frequencies)) {
        frequencies <- x@frequencies
      }
      if (!hasArg(levels)) {
        levels <- intersect(x@levels[[1]], x@levels[[2]])
      }
      if (!hasArg(plotPG)) {
        plotPG <- FALSE
      }
      if (!hasArg(ptw.CIs)) {
        ptw.CIs <- 0.1
      }
      if (!hasArg(type.CIs)) {
        type.CIs <- "naive.sd"
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
                levels.1=levels, levels.2=levels, d1 = 1, d2 = 1) # TODO ... fix
            if (plotPG) {
              PG <- getValues(x@qPG, frequencies = frequencies,
                  levels.1=levels, levels.2=levels, d1 = 1, d2 = 1) # TODO ... fix
            }
            if (hasArg(qsd)) {
              j.min <- round(min(frequencies*2^8/(2*pi)))
              j.max <- round(max(frequencies*2^8/(2*pi)))
              freq.csd <- 2*pi*(j.min:j.max)/2^8
              csd <- getValues(qsd, frequencies = freq.csd,
                  levels.1=levels, levels.2=levels)
            }
            
            #text.headline <- x@weight@descr
            if (ptw.CIs > 0) {
              CI <- getPointwiseCIs(x, frequencies = frequencies,
                  alpha=ptw.CIs, type=type.CIs,
                  levels.1=levels, levels.2=levels, d1 = 1, d2 = 1)
              lowerCIs  <- CI$lowerCIs
              upperCIs  <- CI$upperCIs
              #text.headline <- (paste(text.headline, ", includes ",1-ptw.CIs,"-CI (ptw. of type '",type.CIs,"')",sep=""))
            }
            
            X <- frequencies/(2*pi)
            
            allVals <- array(values[,,,1], dim=c(length(X), K, K))
            if (plotPG)  {
              allVals <- abind(allVals, array(PG[,,,1], dim=c(length(X), K, K)), along=1)
            }
            if (hasArg(qsd)) {
              allVals <- abind(allVals, csd, along=1)
            }
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
                  if (plotPG) {
                    lines(x=X, y=Re(PG[,i1,i2,1]), col=gray(0.5))
                  }
                  if (hasArg(qsd)) {
                    lines(x=freq.csd/(2*pi), y=Re(csd[,i1,i2]), col="red")
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
                  if (plotPG) {
                    lines(x=X, y=Im(PG[,i1,i2,1]), col=gray(0.5))
                  }
                  if (hasArg(qsd)) {
                    lines(x=freq.csd/(2*pi), y=Im(csd[,i1,i2]), col="red")
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
