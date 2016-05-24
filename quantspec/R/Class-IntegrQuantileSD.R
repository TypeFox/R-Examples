#' @include generics.R
#' @include Class-QSpecQuantity.R
#' @include Class-QuantileSD.R
NULL

################################################################################
#' Class for a simulated integrated quantile (i. e., Laplace or copula)
#' density kernel.
#'
#' \code{IntegrQuantileSD} is an S4 class that implements the necessary
#' calculations to determine an integrated version of the quantile spectral
#' density kernel (computed via \code{\link{QuantileSD}}).
#' In particular it can be determined for any model from which a time series
#' of length \code{N} can be sampled via a function call \code{ts(N)}.
#'
#' In the simulation the quantile spectral density is first determined via
#' \code{\link{QuantileSD}}, it's values are recovered using
#' \code{\link{getValues-QuantileSD}} and then cumulated using \code{cumsum}.
#'
#' Note that, all remarks made in the documentation of the super-class
#' \code{\link{QSpecQuantity}} apply.
#'
#' @name   IntegrQuantileSD-class
#' @aliases IntegrQuantileSD
#' @exportClass IntegrQuantileSD
#'
#' @keywords S4-classes
#'
#' @slot qsd a \code{\link{QuantileSD}} from which to begin the computations.
#'
#' @example
#' inst/examples/SmoothedPG-SpecDistr.R
################################################################################

setClass(
    Class = "IntegrQuantileSD",
    representation=representation(
        qsd = "QuantileSD"
    ),
    contains = "QSpecQuantity"
)

setMethod(
    f = "initialize",
    signature = "IntegrQuantileSD",
    definition = function(.Object, qsd) {
      cat("~~~ IntegrQuantileSD: initializator ~~~ \n")

      .Object@qsd <- qsd
      N <- getN(qsd)

      freq <- 2*pi*(0:(N-1))/N
      .Object@frequencies <- freq
      .Object@levels <- getLevels(qsd)

      K1 <- length(getLevels(.Object,1))
      K2 <- length(getLevels(.Object,2))

      .Object@values <- (2*pi/N) * apply(getValues(qsd, frequencies = freq), c(2,3), cumsum)

      return(.Object)
    }
)

################################################################################
#' Get values from a simulated integrated quantile spectral density kernel
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
#' @name getValues-IntegrQuantileSD
#' @aliases getValues,IntegrQuantileSD-method
#'
#' @keywords Access-functions
#'
#' @param object \code{IntegrQuantileSD} of which to get the values
#' @param frequencies a vector of frequencies for which to get the values
#' @param levels.1 the first vector of levels for which to get the values
#' @param levels.2 the second vector of levels for which to get the values
#'
#' @return Returns data from the array \code{values} that's a slot of
#'          \code{object}.
#'
#' @seealso
#' For examples on how to use this function go to \code{\link{IntegrQuantileSD}}.
################################################################################
setMethod(f = "getValues",
    signature = signature("IntegrQuantileSD"),
    definition = function(object,
        frequencies=2*pi*(0:(getN(object@qsd)-1))/getN(object@qsd),
        levels.1=getLevels(object,1),
        levels.2=getLevels(object,2)) {

      # workaround: default values don't seem to work for generic functions?
      if (!hasArg(frequencies)) {
        frequencies <- 2*pi*(0:(getN(object@qsd)-1))/getN(object@qsd)
      }
      if (!hasArg("levels.1")) {
        levels.1 <- getLevels(object,1)
      }
      if (!hasArg("levels.2")) {
        levels.2 <- getLevels(object,2)
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
      r.pos <- closest.pos(oF, f)

      return(object@values[r.pos,c.1.pos,c.2.pos])
    }
)

################################################################################
#' Get associated \code{\link{getQuantileSD}} from an
#' \code{\link{IntegrQuantileSD}}.
#'
#' @name getQuantileSD-IntegrQuantileSD
#' @aliases getQuantileSD,IntegrQuantileSD-method
#'
#' @keywords Access-association-functions
#'
#' @param object \code{IntegrQuantileSD} from which to get the
#'                \code{\link{getQuantileSD}}.
#' @return Returns the \code{\link{getQuantileSD}} object associated.
################################################################################
setMethod(f = "getQuantileSD",
    signature = "IntegrQuantileSD",
    definition = function(object) {
      return(object@qsd)
    }
)

################################################################################
#' Create an instance of the \code{\link{IntegrQuantileSD}} class.
#'
#' @name IntegrQuantileSD-constructor
#' @aliases integrQuantileSD
#' @export
#' 
#' @importFrom stats rnorm
#'
#' @keywords Constructors
#'
#' @param object the number \code{N} of Fourier frequencies to be used;
#'           alternatively a  \code{\link{QuantileSD}} object can be supplied
#'           (then all the other  parameters will be ignored)
#' @param type can be either \code{"Laplace"} or \code{"copula"}; indicates
#'              whether the marginals are to be assumed uniform \eqn{[0,1]}
#'              distributed.
#' @param ts a function that has one argument \code{n} and, each time it is
#'            invoked, returns a new time series from the model for which the
#'            integrated quantile spectral density kernel is to be simulated.
#' @param seed.init an integer serving as an initial seed for the simulations.
#' @param levels.1 A vector of length \code{K1} containing the levels \code{x1}
#'                  at which the \code{QuantileSD} is to be determined.
#' @param levels.2 A vector of length \code{K2} containing the levels \code{x2}
#'                  at which the \code{QuantileSD} is to be determined.
#' @param R an integer that determines the number of independent simulations;
#'           the larger this number the more precise is the result.
#' @param quiet Don't report progress to console when computing the \code{R}
#'               independent quantile periodograms.
#'
#' @return Returns an instance of \code{\link{IntegrQuantileSD}}.
#'
#' @seealso
#' For an example see \code{\link{IntegrQuantileSD}}.
################################################################################
integrQuantileSD <- function(
        object=2^8,
        type=c("copula", "Laplace"),
        ts = rnorm,
        seed.init = 2581,
        levels.1 = 0.5, levels.2=levels.1, R = 1, quiet = FALSE) {

  if (class(object) == "numeric") {
    N <- object
    qsd <- quantileSD(N=N, type=type, ts=ts, seed.init=seed.init,
        levels.1=levels.1, levels.2=levels.2, R=R, quiet=quiet)

  } else if (class(object) == "QuantileSD") {
    qsd <- object
  } else {
    stop("object is neither 'numeric', nor 'QuantileSD'.")
  }


  obj <- new(
      Class = "IntegrQuantileSD",
      qsd = qsd
  )

  return(obj)
}


################################################################################
#' Plot the values of the \code{\link{IntegrQuantileSD}}.
#'
#' Creates a \code{K} x \code{K} plot depicting an integrated quantile spectral
#' density.
#' In each of the subplots either the real part (on and below the diagonal;
#' i. e., \eqn{\tau_1 \leq \tau_2}{tau1 <= tau2}) or the imaginary part
#' (above the diagonal; i. e., \eqn{\tau_1 > \tau_2}{tau1 > tau2}) of
#' \itemize{
#'   \item the integrated quantile spectral density (black line),
#' }
#' for the combination of levels \eqn{\tau_1}{tau1} and \eqn{\tau_2}{tau2}
#' denoted on the left and bottom margin of the plot are displayed.
#'
#' @name plot-IntegrQuantileSD
#' @aliases plot,IntegrQuantileSD,ANY-method
#' @export
#'
#' @param x  The \code{\link{IntegrQuantileSD}} to plot
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
#' @return Plots the simulated integrated quantile spectral density for all
#'          \code{frequencies} and \code{levels} specified.
################################################################################
setMethod(f = "plot",
    signature = signature(x = "IntegrQuantileSD"),
    definition = function(x,
        ratio = 3/2, widthlab = lcm(1), xlab = expression(omega/2*pi), ylab = NULL,
        frequencies=2*pi*(1:(floor(getN(getQuantileSD(x))/2)))/getN(getQuantileSD(x)),
        levels=getLevels(x,1)) {

    def.par <- par(no.readonly = TRUE) # save default, for resetting...

    N <- getN(getQuantileSD(x))

    # workaround: default values don't seem to work for generic functions?
    if (!hasArg(frequencies)) {
      frequencies <- 2*pi*(1:(floor(N/2)))/N
    }
    if (!hasArg(levels)) {
      levels <- getLevels(x,1)
    }
    # end: workaround

		if (length(levels) == 0) {
			stop("There has to be at least one level to plot.")
		}

tryCatch({

    K <- length(levels)
    values <- getValues(x, frequencies = frequencies,
                        levels.1=levels, levels.2=levels)

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
          plot(x=frequencies/(2*pi), y=Re(values[,i1,i2]),
              type="l", xlab="", ylab="")
        } else {
          plot(x=frequencies/(2*pi), y=Im(values[,i1,i2]),
              type="l", xlab="", ylab="")
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
