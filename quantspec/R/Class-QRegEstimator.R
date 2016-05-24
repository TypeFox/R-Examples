#' @include generics.R
#' @include Class-FreqRep.R
#' @include Class-BootPos.R
NULL

################################################################################
#' Class for quantile regression-based estimates in the harmonic linear model.
#'
#' \code{QRegEstimator} is an S4 class that implements the necessary
#' calculations to determine the frequency representation based on the weigthed
#' \eqn{L_1}{L1}-projection of a time series as described in
#' Dette et. al (2015). As a subclass to \code{\link{FreqRep}}
#' it inherits slots and methods defined there.
#'
#' For each frequency \eqn{\omega}{w} from \code{frequencies} and level
#' \eqn{\tau}{tau} from \code{levels} the statistic
#' \deqn{\hat b^{\tau}_n(\omega) := \arg\max_{a \in R, b \in C}
#'       \sum_{t=0}^{n-1}
#'       \rho_{\tau}(Y_t - a - Re(b) \cos(\omega t) - Im(b) \sin(\omega t)),}
#' is determined and stored to the array \code{values}.
#'
#' The solution to the minimization problem is determined using the function
#' \code{\link[quantreg]{rq}} from the \pkg{quantreg} package.
#'
#' All remarks made in the documentation of the super-class
#' \code{\link{FreqRep}} apply.
#'
#' @name   QRegEstimator-class
#' @aliases QRegEstimator
#' @exportClass QRegEstimator
#'
#' @keywords S4-classes
#'
#' @slot method   method used for computing the quantile regression estimates.
#'                 The choice is passed to \code{qr}; see the
#'                 documentation of \code{quantreg} for details.
#' @slot parallel a flag that signalizes that parallelization mechanisms from
#'                 the package \pkg{snowfall} may be used.
#'
#' @references
#' Dette, H., Hallin, M., Kley, T. & Volgushev, S. (2015).
#' Of Copulas, Quantiles, Ranks and Spectra: an \eqn{L_1}{L1}-approach to
#' spectral analysis. \emph{Bernoulli}, \bold{21}(2), 781--831.
#' [cf. \url{http://arxiv.org/abs/1111.7205}]
################################################################################

setClass(
    Class = "QRegEstimator",
    representation=representation(
        method = "character",
        parallel = "logical"
    ),
    contains = "FreqRep"
)

#' @importFrom quantreg rq
#' @importFrom snowfall sfLapply
setMethod(
    f = "initialize",
    signature = "QRegEstimator",
    definition = function(.Object, Y, isRankBased, levels, frequencies, positions.boot, B, method, parallel) {

      .Object@Y <- Y
      .Object@isRankBased <- isRankBased
      .Object@levels <- levels
      .Object@frequencies <- frequencies
      .Object@positions.boot <- positions.boot
      .Object@B <- B
      .Object@method <- method
      .Object@parallel <- parallel

      # Define variables with dimensions
      T <- lenTS(Y)
      K <- length(levels)
      J <- length(frequencies)

      # values[,,,1] contains the non-bootstrapped values
      values <- array(dim=c(J,dim(Y)[2],K,B+1))

      # Convert Y to "pseudo data", if isRankBased == TRUE
      if (isRankBased) {
        data <- apply(Y,2,rank) / T
      } else {
        data <- Y
      }

      qRegSol <- function(X,omega) {
        # Define the harmonic regressors.
        n <- length(X)
        D <- cos(omega*1:n)
        S <- -1*sin(omega*1:n)

        # Then perform the quantile regression.

        suppressWarnings({ # otherwise rq from package quantreg will issue
              # warnings due to non uniqueness of the minimizer.
              if ( abs(omega %% (2*pi)) < .Machine$double.eps^0.5 ) {
                qregSol <- coef(rq(X ~ 1, levels, method=method))
                if (length(levels) > 1) {
                  qregSol <- n * qregSol[1,]
                } else {
                  qregSol <- n * matrix(qregSol[1],ncol=1)
                }
              } else if ( abs((omega-pi) %% (2*pi)) < .Machine$double.eps^0.5 ) {
                qregSol <- coef(rq(X ~ 1 + D, levels, method=method))
                if (length(levels) > 1) {
                  qregSol <- n * qregSol[2,]
                } else {
                  qregSol <- n * matrix(qregSol[2],ncol=1)
                }
              } else {
                qregSol <- coef(rq(X ~ 1 + D + S, levels, method=method))
                pos <- c(2,3)
                iVec <- matrix(c(1, complex(real = 0, imaginary = 1)), nrow = 1)
                if (length(levels)>1) {
                  qregSol <- (n/2) * iVec %*% qregSol[pos,]
                } else {
                  qregSol <- (n/2) * iVec %*% matrix(qregSol[pos],ncol=1)
                }
              }
            })
        return(qregSol)
      }

      for (b in 0:B) {
        for (d in 1:(dim(data)[2])) {
          if (b == 0) {
            qRegSolX <- function(omega){qRegSol(data[,d],omega)}
          } else {
            pos.boot <- getPositions(.Object@positions.boot,B)
            qRegSolX <- function(omega){qRegSol(data[pos.boot,d],omega)}
          }
  
          if (parallel) {
            listVals <- sfLapply(frequencies, qRegSolX)
          } else {
            listVals <- lapply(frequencies, qRegSolX)
          }
          values[,d,,b+1] <- aperm(array(unlist(listVals),dim=c(K,J)), perm=c(2,1))
        }
      }

      #.Object@values <- values
      .Object@values <- array(values, dim=c(J,dim(Y)[2],K,B+1))

      # Return object
      return(.Object)
    }
)

################################################################################
#' Get \code{getParallel} from a \code{\link{QRegEstimator}} object
#'
#' @name getParallel-QRegEstimator
#' @aliases getParallel,QRegEstimator-method
#'
#' @keywords Access-functions
#'
#' @param object \code{QRegEstimator} of which to get the \code{parallel}
#'
#' @return Returns the attribute \code{parallel} that's a slot of \code{object}.
################################################################################
setMethod(f = "getParallel",
    signature = signature("QRegEstimator"),
    definition = function(object) {
      return(object@parallel)
    }
)

################################################################################
#' Create an instance of the \code{QRegEstimator} class.
#'
#' The parameter \code{type.boot} can be set to choose a block bootstrapping
#' procedure. If \code{"none"} is chosen, a moving blocks bootstrap with
#' \code{l=length(Y)} and  \code{N=length(Y)} would be done. Note that in that
#' case one would also chose \code{B=0} which means that \code{getPositions}
#' would never be called. If \code{B>0} then each bootstrap replication would
#' be the undisturbed time series.
#'
#' @name QRegEstimator-constructor
#' @aliases qRegEstimator
#' @export
#' @importFrom snowfall sfInit sfLibrary sfExportAll sfStop
#'
#' @keywords Constructors
#'
#' @param Y A \code{vector} of real numbers containing the time series from
#'           which to determine the quantile periodogram or a \code{ts} object
#'           or a \code{zoo} object.
#' @param isRankBased If true the time series is first transformed to pseudo
#'                    data [cf. \code{\link{FreqRep}}].
#' @param levels A vector of length \code{K} containing the levels \code{x}
#'                  at which the \code{QRegEstimator} is to be determined.
#' @param frequencies A vector containing frequencies at which to determine the
#'                     \code{QRegEstimator}.
#' @param B number of bootstrap replications
#' @param l (expected) length of blocks
#' @param type.boot A flag to choose a method for the block bootstrap; currently
#'                   two options are implemented: \code{"none"} and \code{"mbb"}
#'                   which means to do a moving blocks  bootstrap with \code{B}
#'                   and \code{l} as specified.
#' @param method  method used for computing the quantile regression estimates.
#'                 The choice is passed to \code{qr}; see the
#'                 documentation of \code{quantreg} for details.
#' @param parallel a flag to allow performing parallel computations.
#'
#' @return Returns an instance of \code{QRegEstimator}.
#'
#' @example
#' inst/examples/QRegEstimator-parallel.R
################################################################################
qRegEstimator <- function( Y,
    frequencies=2*pi/lenTS(Y) * 0:(lenTS(Y)-1),
    levels = 0.5,
    isRankBased=TRUE,
    B = 0,
    l = 0,
    type.boot = c("none", "mbb"),
    method = c("br", "fn", "pfn", "fnc", "lasso", "scad"),
    parallel = FALSE) {

  # Verify if all parameters are valid
  Y <- timeSeriesValidator(Y)

  if (!(is.vector(frequencies)  && is.numeric(frequencies))) {
    stop("'frequencies' needs to be specified as a vector of real numbers")
  }

  if (!(is.vector(levels) && is.numeric(levels))) {
    stop("'levels' needs to be specified as a vector of real numbers")
  }

  if (isRankBased && !(prod(levels >= 0) && prod(levels <=1))) {
    stop("'levels' need to be from [0,1] when isRankBased==TRUE")
  }

  # Check validity of frequencies
  frequencies <- frequenciesValidator(frequencies, lenTS(Y))

  type.boot <- match.arg(type.boot, c("none","mbb"))[1]
  switch(type.boot,
      "none" = {
        bootPos <- movingBlocks(lenTS(Y),lenTS(Y))},
      "mbb" = {
        bootPos <- movingBlocks(l,lenTS(Y))}
  )

  method <- match.arg(method, c("br", "fn", "pfn", "fnc", "lasso", "scad"))[1]

  freqRep <- new(
      Class = "QRegEstimator",
      Y = Y,
      isRankBased = isRankBased,
      levels = sort(levels),
      B = B,
      positions.boot = bootPos,
      frequencies = frequencies,
      method = method,
      parallel = parallel
  )

  return(freqRep)
}
