#' @include generics.R
#' @include Class-QSpecQuantity.R
NULL

################################################################################
#' Class for a simulated quantile (i. e., Laplace or copula)
#' density kernel.
#'
#' \code{QuantileSD} is an S4 class that implements the necessary
#' calculations to determine a numeric approximation to the quantile spectral
#' density kernel of a model from which a time series of length \code{N} can be
#' sampled via a function call \code{ts(N)}.
#'
#' In the simulation a number of \code{R} independent quantile periodograms
#' based on the clipped time series are simulated. If \code{type=="copula"},
#' then the rank-based version is used. The sum and the sum of the squared
#' absolute value is stored to the slots \code{sumPG} and \code{sumSqPG}.
#' After the simulation is completed the mean and it's standard error (of the
#' simulated quantile periodograms) are determined and stored to \code{meanPG}
#' and \code{stdError}. Finally, the (copula) spectral density kernel is
#' determined by smoothing real and imaginary part of \code{meanPG} seperately
#' for each combination of levels using \code{\link[stats]{smooth.spline}}.
#'
#' Note that, all remarks made in the documentation of the super-class
#' \code{\link{QSpecQuantity}} apply.
#'
#' @name   QuantileSD-class
#' @aliases QuantileSD
#' @exportClass QuantileSD
#'
#' @keywords S4-classes
#'
#' @slot N a \code{numeric} specifying the number of equaly spaced
#'           Fourier frequencies from
#'           \eqn{[0,2\pi)}{[0,2pi)} for which the (copula) spectral density
#'           will be simulated; note that due to the simulation mechanism a
#'           larger number will also yield a better approximation.
#' @slot R the number of independent repetitions performed; note that due to
#'          the simulation mechanism a larger number will also yield a better
#'          approximation; can be enlarged
#'          using \code{\link{increasePrecision-QuantileSD}}.
#' @slot type can be either \code{Laplace} or \code{copula}; indicates whether
#'             the marginals are to be assumed uniform \eqn{[0,1]} distributed.
#' @slot ts a \code{function} that allows to draw independent samples
#'           \eqn{Y_0, \ldots, Y_{n-1}} from the process for which the (copula)
#'           spectral density kernel is to be simulated
#' @slot seed.last used internally to store the state of the pseudo random number
#'                  generator, so the precision can be increased by generating
#'                  more pseudo random numbers that are independent from the ones
#'                  previously used.
#' @slot sumPG an \code{array} used to store the sum of the simulated quantile
#'              periodograms
#' @slot sumSqPG an \code{array} used to store the sum of the squared absolute
#'                values of the simulated quantile periodograms
#' @slot meanPG an \code{array} used to store the mean of the simulated quantile
#'               periodograms
#' @slot stdError an \code{array} used to store the estimated standard error of the mean
#'                  of the simulated quantile periodograms
#'
#' @seealso
#' Examples for implementations of functions \code{ts} can be found at:
#' \code{\link{ts-models}}.
#'
#' @references
#' Dette, H., Hallin, M., Kley, T. & Volgushev, S. (2015).
#' Of Copulas, Quantiles, Ranks and Spectra: an \eqn{L_1}{L1}-approach to
#' spectral analysis. \emph{Bernoulli}, \bold{21}(2), 781--831.
#' [cf. \url{http://arxiv.org/abs/1111.7205}]
#'
#' Kley, T., Volgushev, S., Dette, H. & Hallin, M. (2016).
#' Quantile Spectral Processes: Asymptotic Analysis and Inference.
#' \emph{Bernoulli}, \bold{22}(3), 1770--1807.
#' [cf. \url{http://arxiv.org/abs/1401.8104}]
#' 
#' Barunik, J. & Kley, T. (2015).
#' Quantile Cross-Spectral Measures of Dependence between Economic Variables.
#' [preprint available from the authors]
#'
#' @example
#' inst/examples/QuantileSD.R
################################################################################
setClass(
    Class = "QuantileSD",
    representation=representation(
        N = "numeric",
        R = "numeric",
        type = "character",
        ts = "function",
        seed.last = "list",
        sumPG = "array",
        sumSqPG = "array",
        meanPG = "array",
        stdError = "array"
    ),
    contains = "QSpecQuantity"
)

################################################################################
#' Increase the precision of a \code{QuantileSD}
#'
#' The precision is increased by generating an additional \code{R}
#' \code{\link{QuantilePG}} objects (independent of the previous ones) and
#' then including them in the average.
#'
#' @name increasePrecision-QuantileSD
#' @aliases increasePrecision,QuantileSD-method
#' 
#' @importFrom stats smooth.spline
#'
#' @param object The \code{\link{QuantileSD}} of which to increase the precision.
#' @param R value of which to enlarge R
#' @param quiet Don't report progress to console when computing the \code{R}
#'               independent quantile periodograms.
#'
#' @return Returns an \code{\link{QuantileSD}} object determined from
#'          \code{oldR + R} independent repetitions.
#'
#' @examples
#' # First simulate a copula spectral density from R=20 independent runs.
#' csd <- quantileSD(N=2^9, ts=ts1, levels.1=c(0.25,0.5), type="copula", R=20)
#'
#' # Check out the result:
#' getR(csd)
#' plot(csd)
#'
#' # Now increase the number of independent simulation runs to 50.
#' csd <- increasePrecision(csd, R=30)
#'
#' # Check out the (more precise) result:
#' getR(csd)
#' plot(csd)
################################################################################
setMethod(f = "increasePrecision",
    signature = signature("QuantileSD"),
    definition = function(object, R=1, quiet = FALSE) {

      if (!hasArg(R)) {
        R <- 1
      }

      ts <- object@ts
      N <- object@N
      L <- floor(N/2)
      levels.1 <- object@levels[[1]]
      levels.2 <- object@levels[[2]]

      sumPG <- object@sumPG
      sumSqPG <- object@sumSqPG
      
      J <- dim(sumPG)[1]
      P <- dim(sumPG)[2]

      type <- object@type

      switch(type,
          "copula" = {
            rb <- TRUE},
          "Laplace" = {
            rb <- FALSE}
      )


      save_rng <- function(savefile=tempfile()) {
        if (exists(".Random.seed"))  {
          oldseed <- get(".Random.seed", .GlobalEnv)
        } else stop("don't know how to save before set.seed() or r*** call")
        oldRNGkind <- RNGkind()

        return(list(seed=oldseed,RNGkind=oldRNGkind))
      }

      restore_rng <- function(savedseed) {
        do.call("RNGkind",as.list(savedseed$RNGkind))  ## must be first!
        assign(".Random.seed", savedseed$seed, .GlobalEnv)
      }

      # Begin of Program

      restore_rng(object@seed.last)

      K1 <- length(levels.1)
      K2 <- length(levels.2)

      progress <- 0
      beginTime <- Sys.time()
      freq <- object@frequencies
      freq[which(freq == 0)] <- freq[which(freq == min(freq[freq>0]))]
      for (r in 1:R) {
        if (!quiet && (round(r/R*100) != progress)) {
          progress <- round(r/R*100)
          timeElapsed <- difftime(time1=Sys.time(), time2=beginTime)
          cat(paste(progress,"% done; ", format(timeElapsed, format="%H:%M:%S"), "elapsed \n"))
        }

        Y <- ts(N)

        A <- getValues(quantilePG(Y, type="clipped", isRankBased=rb,
                levels.1 = levels.1, levels.2 = levels.2),
            frequencies = freq)
        if (length(dim(A)) == 4) {
          A <- A[,,,1]
        } else {
          A <- A[,,,,,1]
        }
        A <- array(A, dim = c(J,P,K1,P,K2))
        sumPG <- sumPG + A 
        sumSqPG <- sumSqPG + abs(A)^2
      }

      object@R <- object@R + R
      object@stdError <- Re(1/(object@R-1) * (sumSqPG - 1/(object@R)*abs(sumPG)^2)) / sqrt(object@R)

      object@sumPG <- sumPG
      object@sumSqPG <- sumSqPG
      meanPG <- sumPG / (object@R)
      object@meanPG <- meanPG


      values <- array(0, dim=c(L+1,P,K1,P,K2))
      for (p1 in 1:P) {
        for (p2 in 1:P) {
          for (j1 in 1:K1) {
            for (j2 in 1:K2) {
              Rp <- smooth.spline(x=0:L, y=Re(meanPG[,p1,j1,p2,j2]))$y
              Ip <- smooth.spline(x=0:L, y=Im(meanPG[,p1,j1,p2,j2]))$y
              values[,p1,j1,p2,j2] <- complex(real=Rp, imaginary=Ip)
            }
          }
        }
      }
      object@values <- values
      object@seed.last <- save_rng()

      # End of main Program

      return(object)

    }
)

setMethod(
    f = "initialize",
    signature = "QuantileSD",
    definition = function(.Object, type, N, levels, ts, R, seed.last, quiet) {

      L <- floor(N/2)
      freq <- 2*pi*(0:L)/N

      .Object@N <- N
      .Object@frequencies <- freq
      .Object@levels <- levels
      .Object@ts <- ts
      .Object@R <- 0
      .Object@type <- type

      K1 <- length(getLevels(.Object,1))
      K2 <- length(getLevels(.Object,2))
      
      #TODO: check the next two lines again!
      Y_test <- timeSeriesValidator(ts(2))
      P <- dim(Y_test)[2]

      .Object@sumPG <- array(0, dim=c(L+1,P,K1,P,K2))
      .Object@sumSqPG <- array(0, dim=c(L+1,P,K1,P,K2))

      .Object@seed.last <- seed.last

      .Object <- increasePrecision(.Object, R, quiet)

      return(.Object)
    }
)


################################################################################
#' Get values from a quantile spectral density kernel
#'
#' If none of the optional parameters is specified then the values are returned
#' for all Fourier frequencies in \eqn{[0,2\pi)}{[0,2pi)} (base given by slot
#' \code{N}) and all levels available. The frequencies and levels can be freely
#' specified. The returned array then has, at position \code{(j,k1,k2,b)},
#' the value corresponding to the \code{frequencies[j]},
#' \code{levels.1[k1]} and \code{levels.2[k2]} that are closest to the
#' \code{frequencies}, \code{levels.1} and \code{levels.2}
#' available in \code{object}; \code{\link{closest.pos}} is used to determine
#' what closest to means.
#'
#' @name getValues-QuantileSD
#' @aliases getValues,QuantileSD-method
#'
#' @keywords Access-functions
#'
#' @param object \code{QuantileSD} of which to get the values
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
#' For examples on how to use this function go to \code{\link{QuantileSD}}.
################################################################################
setMethod(f = "getValues",
    signature = signature("QuantileSD"),
    definition = function(object,
        frequencies=2*pi*(0:(object@N-1))/object@N,
        levels.1=getLevels(object,1),
        levels.2=getLevels(object,2),
        d1 = 1:(dim(object@values)[2]),
        d2 = 1:(dim(object@values)[4])) {

      # workaround: default values don't seem to work for generic functions?
      if (!hasArg(frequencies)) {
        frequencies <- 2*pi*(0:(object@N-1))/object@N
      }
      if (!hasArg("levels.1")) {
        levels.1 <- getLevels(object,1)
      }
      if (!hasArg("levels.2")) {
        levels.2 <- getLevels(object,2)
      }
      if (!hasArg("d1")) {
        d1 <- 1:(dim(object@values)[2])
      }
      if (!hasArg("d2")) {
        d2 <- 1:(dim(object@values)[4])
      }
      # end: workaround

      ##############################
      ## (Similar) Code also in Class-FreqRep!!!
      ## (Similar) Code also in Class-SmoothedPG!!!
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
      c.1.pos <- closest.pos(getLevels(object,1),levels.1)
      c.2.pos <- closest.pos(getLevels(object,2),levels.2)

      if (!subsetequal.approx(levels.1, getLevels(object,1))) {
        warning("Not all 'values' for 'levels.1' requested were available. 'values' for the next available level are returned.")
      }

      if (!subsetequal.approx(levels.2, getLevels(object,2))) {
        warning("Not all 'values' for 'levels.2' requested were available. 'values' for the next available level are returned.")
      }

      # Select rows
      r1.pos <- closest.pos(oF, f[f <= pi])
      r2.pos <- closest.pos(-1*(2*pi-oF),-1*f[f > pi])

      J <- length(frequencies)
      K1 <- length(levels.1)
      K2 <- length(levels.2)
      D1 <- length(d1)
      D2 <- length(d2)
      res <- array(dim=c(J, D1, K1, D2, K2))

      if (length(r1.pos) > 0) {
        res[which(f <= pi),,,,] <- object@values[r1.pos,d1,c.1.pos,d2,c.2.pos]
      }
      if (length(r2.pos) > 0) {
        res[which(f > pi),,,,] <- Conj(object@values[r2.pos,d1,c.1.pos,d2,c.2.pos])
      }
      

      if (D1 == 1 && D2 == 1) {
        final.dim.res <- c(J, K1, K2)
      } else {
        final.dim.res <- c(J, D1, K1, D2, K2)
      }
      
      res <- array(res, dim=final.dim.res)

      return(res)
    }
)


################################################################################
#' Compute quantile coherency from a quantile spectral density kernel
#'
#' Returns quantile coherency defined as
#' \deqn{\frac{f^{j_1, j_2}(\omega; \tau_1, \tau_2)}{(f^{j_1, j_1}(\omega; \tau_1, \tau_1) f^{j_2, j_2}(\omega; \tau_2, \tau_2))^{1/2}}}
#' where \eqn{f^{j_1, j_2}(\omega; \tau_1, \tau_2)} is the quantile spectral density.
#' 
#' For the mechanism of selecting frequencies, dimensions and/or levels see,
#' for example, \code{\link{getValues-QuantileSD}}.
#' 
#' @name getCoherency-QuantileSD
#' @aliases getCoherency,QuantileSD-method
#'
#' @keywords Access-functions
#'
#' @param object \code{QuantileSD} of which to get the values
#' @param frequencies a vector of frequencies for which to get the values
#' @param levels.1 the first vector of levels for which to get the values
#' @param levels.2 the second vector of levels for which to get the values
#' @param d1 optional parameter that determine for which j1 to return the
#' 					 data; may be a vector of elements 1, ..., D
#' @param d2 same as d1, but for j2
#'
#' @return Returns data from the coherency as defined in the details.
#'
#' @seealso
#' For examples on how to use this function go to \code{\link{QuantileSD}}.
################################################################################
setMethod(f = "getCoherency",
    signature = signature("QuantileSD"),
    definition = function(object,
        frequencies=2*pi*(0:(object@N-1))/object@N,
        levels.1=getLevels(object,1),
        levels.2=getLevels(object,2),
        d1 = 1:(dim(object@values)[2]),
        d2 = 1:(dim(object@values)[4])) {
      
      # workaround: default values don't seem to work for generic functions?
      if (!hasArg(frequencies)) {
        frequencies <- 2*pi*(0:(object@N-1))/object@N
      }
      if (!hasArg("levels.1")) {
        levels.1 <- getLevels(object,1)
      }
      if (!hasArg("levels.2")) {
        levels.2 <- getLevels(object,2)
      }
      if (!hasArg("d1")) {
        d1 <- 1:(dim(object@values)[2])
      }
      if (!hasArg("d2")) {
        d2 <- 1:(dim(object@values)[4])
      }
      # end: workaround
      
#      ##############################
#      ## (Similar) Code also in Class-FreqRep!!!
#      ## (Similar) Code also in Class-SmoothedPG!!!
#      ##############################
#      
#      # Transform all frequencies to [0,2pi)
#      frequencies <- frequencies %% (2*pi)
#      
#      # Create an aux vector with all available frequencies
#      oF <- object@frequencies
#      f <- frequencies
#      
#      # returns TRUE if x c y
#      subsetequal.approx <- function(x,y) {
#        X <- round(x, .Machine$double.exponent-2)
#        Y <- round(y, .Machine$double.exponent-2)
#        return(setequal(X,intersect(X,Y)))
#      }
#      
#      C1 <- subsetequal.approx(f[f <= pi], oF)
#      C2 <- subsetequal.approx(f[f > pi], 2*pi - oF[which(oF != 0 & oF != pi)])
#      
#      if (!(C1 & C2)) {
#        warning("Not all 'values' for 'frequencies' requested were available. 'values' for the next available Fourier frequencies are returned.")
#      }
#      
#      # Select columns
#      c.1.pos <- closest.pos(getLevels(object,1),levels.1)
#      c.2.pos <- closest.pos(getLevels(object,2),levels.2)
#      
#      if (!subsetequal.approx(levels.1, getLevels(object,1))) {
#        warning("Not all 'values' for 'levels.1' requested were available. 'values' for the next available level are returned.")
#      }
#      
#      if (!subsetequal.approx(levels.2, getLevels(object,2))) {
#        warning("Not all 'values' for 'levels.2' requested were available. 'values' for the next available level are returned.")
#      }
#      
#      # Select rows
#      r1.pos <- closest.pos(oF, f[f <= pi])
#      r2.pos <- closest.pos(-1*(2*pi-oF),-1*f[f > pi])
      
      J <- length(frequencies)
      K1 <- length(levels.1)
      K2 <- length(levels.2)
      D1 <- length(d1)
      D2 <- length(d2)
      res <- array(dim=c(J, D1, K1, D2, K2))
      
#      if (length(r1.pos) > 0) {
#        res[which(f <= pi),,,,] <- object@values[r1.pos,,c.1.pos,,c.2.pos]
#      }
#      if (length(r2.pos) > 0) {
#        res[which(f > pi),,,,] <- Conj(object@values[r2.pos,,c.1.pos,,c.2.pos])
#      }
#

    #if (class(object@weight) != "KernelWeight") {
    #  error("Coherency can only be determined if weight is of type KernelWeight.")
    #}
    d <- union(d1,d2)
    V <- getValues(object, d1 = d, d2 = d, frequencies = frequencies)
    V <- array(V, dim = c(dim(V), 1))
    
    d1.pos <- closest.pos(d, d1)
    d2.pos <- closest.pos(d, d2)
    
    res <- .computeCoherency(V, d1.pos, d2.pos)

      if (D1 == 1 && D2 == 1) {
        final.dim.res <- c(J, K1, K2)
      } else {
        final.dim.res <- c(J, D1, K1, D2, K2)
      }

      res <- array(res, dim=final.dim.res)
      
      return(res)
    }
)

################################################################################
#' Get \code{N} from a quantile spectral density kernel
#'
#' @name getN-QuantileSD
#' @aliases getN,QuantileSD-method
#'
#' @keywords Access-functions
#'
#' @param object \code{QuantileSD} of which to get the \code{N}
#' @return Returns the attribute \code{N} that's a slot of \code{object}.
################################################################################
setMethod(f = "getN",
    signature = signature("QuantileSD"),
    definition = function(object) {
      return(object@N)
    }
)

################################################################################
#' Get \code{R} from a quantile spectral density kernel
#'
#' @name getR-QuantileSD
#' @aliases getR,QuantileSD-method
#'
#' @keywords Access-functions
#'
#' @param object \code{QuantileSD} of which to get the \code{R}
#' @return Returns the attribute \code{R} that's a slot of \code{object}.
################################################################################
setMethod(f = "getR",
    signature = signature("QuantileSD"),
    definition = function(object) {
      return(object@R)
    }
)

################################################################################
#' Get \code{type} from a quantile spectral density kernel
#'
#' @name getType-QuantileSD
#' @aliases getType,QuantileSD-method
#'
#' @keywords Access-functions
#'
#' @param object \code{QuantileSD} of which to get the \code{type}
#' @return Returns the attribute \code{type} that's a slot of
#'         \code{object}.
################################################################################
setMethod(f = "getType",
    signature = signature("QuantileSD"),
    definition = function(object) {
      return(object@type)
    }
)

################################################################################
#' Get \code{ts} from a quantile spectral density kernel
#'
#' @name getTs-QuantileSD
#' @aliases getTs,QuantileSD-method
#'
#' @keywords Access-functions
#'
#' @param object \code{QuantileSD} of which to get the \code{ts}
#' @return Returns the attribute \code{ts} that's a slot of \code{object}.
################################################################################
setMethod(f = "getTs",
    signature = signature("QuantileSD"),
    definition = function(object) {
      return(object@ts)
    }
)

################################################################################
#' Get \code{meanPG} from a quantile spectral density kernel
#'
#' The selection mechanism for frequencies and levels operates in the same way
#' as described in \code{\link{getValues-QuantileSD}}. The format of the
#' output is also described there.
#'
#' @name getMeanPG-QuantileSD
#' @aliases getMeanPG,QuantileSD-method
#'
#' @keywords Access-functions
#'
#' @param object \code{QuantileSD} of which to get the \code{meanPG}
#' @param frequencies a vector of frequencies for which to get the \code{meanPG}
#' @param levels.1 the first vector of levels for which to get the \code{meanPG}
#' @param levels.2 the second vector of levels for which to get the \code{meanPG}
#' @param d1 optional parameter that determine for which j1 to return the
#' 					 \code{meanPG}; may be a vector of elements 1, ..., D
#' @param d2 same as d1, but for j2
#'
#' @return Returns the array \code{meanPG} that's a slot of \code{object}.
################################################################################
setMethod(f = "getMeanPG",
    signature = signature("QuantileSD"),
    definition = function(object,
        frequencies=2*pi*(0:(getN(object)-1))/getN(object),
        levels.1=getLevels(object,1),
        levels.2=getLevels(object,2),
        d1 = 1:(dim(object@values)[2]),
        d2 = 1:(dim(object@values)[4])) {

      # workaround: default values don't seem to work for generic functions?
      if (!hasArg(frequencies)) {
        frequencies <- 2*pi*(0:(object@N-1))/object@N
      }
      if (!hasArg("levels.1")) {
        levels.1 <- getLevels(object,1)
      }
      if (!hasArg("levels.2")) {
        levels.2 <- getLevels(object,2)
      }
      if (!hasArg("d1")) {
        d1 <- 1:(dim(object@values)[2])
      }
      if (!hasArg("d2")) {
        d2 <- 1:(dim(object@values)[4])
      }
      # end: workaround

      pfreq <- closest.pos(object@frequencies, frequencies)
      plevels1 <- closest.pos(getLevels(object,1), levels.1)
      plevels2 <- closest.pos(getLevels(object,2), levels.2)
      
      J <- length(pfreq)
      K1 <- length(plevels1)
      K2 <- length(plevels2)
      D1 <- length(d1)
      D2 <- length(d2)
      
      res <- object@meanPG[pfreq, d1, plevels1, d2, plevels2, drop=F]
      
      final.dim.res <- c(J)
      if (D1 > 1) {
        final.dim.res <- c(final.dim.res,D1)
      }
      final.dim.res <- c(final.dim.res,K1)
      if (D2 > 1) {
        final.dim.res <- c(final.dim.res,D2)
      }
      final.dim.res <- c(final.dim.res,K2)
      
      res <- array(res, dim=final.dim.res)
      
      return(res)
    }
)

################################################################################
#' Get \code{stdError} from a quantile spectral density kernel
#'
#' The selection mechanism for frequencies and levels operates in the same way
#' as described in \code{\link{getValues-QuantileSD}}. The format of the
#' output is also described there.
#'
#' @name getStdError-QuantileSD
#' @aliases getStdError,QuantileSD-method
#'
#' @keywords Access-functions
#'
#' @param object \code{QuantileSD} of which to get the \code{stdError}
#' @param frequencies a vector of frequencies for which to get the \code{stdError}
#' @param levels.1 the first vector of levels for which to get the \code{stdError}
#' @param levels.2 the second vector of levels for which to get the \code{stdError}
#' @param d1 optional parameter that determine for which j1 to return the
#' 					 \code{stdError}; may be a vector of elements 1, ..., D
#' @param d2 same as d1, but for j2
#'
#' @return Returns the array \code{stdError} that's a slot of \code{object}.
################################################################################
setMethod(f = "getStdError",
    signature = signature("QuantileSD"),
    definition = function(object,
        frequencies=2*pi*(0:(object@N-1))/object@N,
        levels.1=getLevels(object,1),
        levels.2=getLevels(object,2),
        d1 = 1:(dim(object@values)[2]),
        d2 = 1:(dim(object@values)[4])) {

      # workaround: default values don't seem to work for generic functions?
      if (!hasArg(frequencies)) {
        frequencies <- 2*pi*(0:(object@N-1))/object@N
      }
      if (!hasArg("levels.1")) {
        levels.1 <- getLevels(object,1)
      }
      if (!hasArg("levels.2")) {
        levels.2 <- getLevels(object,2)
      }
      if (!hasArg("d1")) {
        d1 <- 1:(dim(object@values)[2])
      }
      if (!hasArg("d2")) {
        d2 <- 1:(dim(object@values)[4])
      }
      # end: workaround

      pfreq <- closest.pos(object@frequencies, frequencies)
      plevels1 <- closest.pos(getLevels(object,1), levels.1)
      plevels2 <- closest.pos(getLevels(object,2), levels.2)
      
      J <- length(pfreq)
      K1 <- length(plevels1)
      K2 <- length(plevels2)
      D1 <- length(d1)
      D2 <- length(d2)
      
      res <- object@stdError[pfreq, d1, plevels1, d2, plevels2, drop=F]
      
      final.dim.res <- c(J)
      if (D1 > 1) {
        final.dim.res <- c(final.dim.res,D1)
      }
      final.dim.res <- c(final.dim.res,K1)
      if (D2 > 1) {
        final.dim.res <- c(final.dim.res,D2)
      }
      final.dim.res <- c(final.dim.res,K2)
      
      res <- array(res, dim=final.dim.res)
      
      return(res)

    }
)

################################################################################
#' Get associated \code{\link{QuantilePG}} from a \code{\link{QuantileSD}}.
#'
#' @name getQuantilePG-QuantileSD
#' @aliases getQuantilePG,QuantileSD-method
#'
#' @keywords Access-association-functions
#'
#' @param object \code{QuantileSD} from which to get the
#'                \code{\link{QuantilePG}}.
#'
#' @return Returns the \code{\link{QuantilePG}} object associated.
################################################################################
setMethod(f = "getQuantilePG",
    signature = "QuantileSD",
    definition = function(object) {
      return(object@qPG)
    }
)

################################################################################
#' Create an instance of the \code{\link{QuantileSD}} class.
#'
#' @name QuantileSD-constructor
#' @aliases quantileSD
#' @export
#' 
#' @importFrom stats rnorm
#' @importFrom stats runif
#'
#' @keywords Constructors
#'
#' @param N the number of Fourier frequencies to be used.
#' @param type can be either \code{Laplace} or \code{copula}; indicates whether
#'             the marginals are to be assumed uniform \eqn{[0,1]} distributed.
#' @param ts a function that has one argument \code{n} and, each time it is
#'            invoked, returns a new time series from the model for which the
#'            copula spectral density kernel is to be simulated.
#' @param seed.init an integer serving as an initial seed for the simulations.
#' @param levels.1 A vector of length \code{K1} containing the levels \code{x1}
#'                  at which the \code{QuantileSD} is to be determined.
#' @param levels.2 A vector of length \code{K2} containing the levels \code{x2}
#'                  at which the \code{QuantileSD} is to be determined.
#' @param R an integer that determines the number of independent simulations;
#'           the larger this number the more precise is the result.
#' @param quiet Dont't report progress to console when computing the \code{R}
#'               independent quantile periodograms.
#'
#' @return Returns an instance of \code{QuantileSD}.
#'
#' @seealso
#' For examples see \code{\link{QuantileSD}}.
################################################################################
quantileSD <- function(N=2^8,
    type = c("copula","Laplace"),
    ts = rnorm,
    seed.init = runif(1),
    levels.1, levels.2=levels.1, R = 1, quiet = FALSE) {

  # helper function
  save_rng <- function(savefile=tempfile()) {
    if (exists(".Random.seed"))  {
      oldseed <- get(".Random.seed", .GlobalEnv)
    } else stop("don't know how to save before set.seed() or r*** call")
    oldRNGkind <- RNGkind()

    return(list(seed=oldseed,RNGkind=oldRNGkind))
  }

  set.seed(seed.init)

  type <- match.arg(type)[1]

  obj <- new(
      Class = "QuantileSD",
      type = type,
      N = N,
      levels = list(levels.1, levels.2),
      ts = ts,
      R = R,
      seed.last = save_rng(),
      quiet = quiet
  )

  return(obj)
}

################################################################################
#' Plot the values of the \code{\link{QuantileSD}}.
#'
#' Creates a \code{K} x \code{K} plot depicting a quantile spectral density.
#' In each of the subplots either the real part (on and below the diagonal;
#' i. e., \eqn{\tau_1 \leq \tau_2}{tau1 <= tau2}) or the imaginary parts
#' (above the diagonal; i. e., \eqn{\tau_1 > \tau_2}{tau1 > tau2}) of
#' \itemize{
#'   \item the quantile spectral density (red line),
#'   \item the means of the quantile periodograms used in the simulation
#'         (black line),
#' }
#' for the combination of levels \eqn{\tau_1}{tau1} and \eqn{\tau_2}{tau2}
#' denoted on the left and bottom margin of the plot are displayed.
#' 
#' Currently, only the plot for the first component is shown.
#'
#' @name plot-QuantileSD
#' @aliases plot,QuantileSD,ANY-method
#' @export
#'
#' @param x  The \code{\link{QuantileSD}} to plot
#' @param ratio quotient of width over height of the subplots; use this
#'               parameter to produce landscape or portrait shaped plots.
#' @param widthlab width for the labels (left and bottom); default is
#' 							   \code{lcm(1)}, cf. \code{\link[graphics]{layout}}.
#' @param xlab label that will be shown on the bottom of the plots; can be
#' 						 an expression (for formulas), characters or \code{NULL} to
#' 						 force omission (to save space).
#' @param ylab label that will be shown on the left side of the plots;
#' 						 can be an expression (for formulas), characters or
#' 						 \code{NULL} to force omission (to save space).
#' @param frequencies a set of frequencies for which the values are to be
#'                    plotted.
#' @param levels a set of levels for which the values are to be plotted.
#'
#' @return Plots the simulated quantile spectral density for all
#'          \code{frequencies} and \code{levels} specified.
################################################################################

setMethod(f = "plot",
    signature = signature(x = "QuantileSD"),
    definition = function(x,
        ratio = 3/2, widthlab = lcm(1), xlab = expression(omega/2*pi), ylab = NULL,
        frequencies=2*pi*(1:(floor(x@N/2)))/x@N,
        levels=getLevels(x,1)) {

    def.par <- par(no.readonly = TRUE) # save default, for resetting...

    # workaround: default values don't seem to work for generic functions?
    if (!hasArg(frequencies)) {
      frequencies <- 2*pi*(1:(floor(x@N/2)))/x@N
    }
    if (!hasArg(levels)) {
      levels <- getLevels(x,1)
    }
    # end: workaround

    if (length(levels) == 0) {
      stop("There has to be at least one level to plot.")
    }

tryCatch({

    N <- x@N
    K <- length(levels)
    values <- getValues(x, frequencies = frequencies,
                        levels.1=levels, levels.2=levels, d1=1, d2=1)
    meanPG <- getMeanPG(x, frequencies = frequencies,
        levels.1=levels, levels.2=levels, d1=1, d2=1)

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
    # END TEST

    for (i1 in 1:K) {
      for (i2 in 1:K) {
        if (i2 >= i1) {
          plot(x=frequencies/(2*pi), y=Re(meanPG[,i1,i2]),
              type="l", xlab="", ylab="")
          lines(x=frequencies/(2*pi), y=Re(values[,i1,i2]), col="red")
        } else {
          plot(x=frequencies/(2*pi), y=Im(meanPG[,i1,i2]),
              type="l", xlab="", ylab="")
          lines(x=frequencies/(2*pi), y=Im(values[,i1,i2]), col="red")
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
