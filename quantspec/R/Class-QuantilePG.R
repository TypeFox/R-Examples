#' @include generics.R
#' @include Class-QSpecQuantity.R
#' @include Class-FreqRep.R
NULL

################################################################################
#' Class for a quantile (i. e., Laplace or copula) periodogram.
#'
#' \code{QuantilePG} is an S4 class that implements the necessary
#' calculations to determine one of the periodogram-like statistics defined in
#' Dette et. al (2015) and Kley et. al (2016).
#'
#' Performs all the calculations to determine a quantile periodogram from a
#' \code{FreqRep} object upon initizalization (and on request
#' stores the values for faster access).
#' The two methods available for the estimation are the ones implemented as
#' subclasses of \code{\link{FreqRep}}:
#' \itemize{
#'     \item the Fourier transformation of the clipped time series
#'           \eqn{(\{I\{Y_t \leq q\})}{(I{Y_t <= q})} [cf. \code{\link{ClippedFT}}], or
#'     \item the weighted \eqn{L_1}{L1}-projection of \eqn{(Y_t)} onto an harmonic
#'           basis [cf. \code{\link{QRegEstimator}}].
#' }
#'
#' All remarks made in the documentation of the super-class
#' \code{\link{QSpecQuantity}} apply.
#'
#' @name   QuantilePG-class
#' @aliases QuantilePG
#' @exportClass QuantilePG
#'
#' @keywords S4-classes
#'
#' @slot freqRep a \code{\link{FreqRep}} object where the quantile periodogram will
#'       be based on.
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
#' @example
#' inst/examples/QuantilePG.R
################################################################################
setClass(
    Class = "QuantilePG",
    representation=representation(
        freqRep = "FreqRep"
    ),
    contains = "QSpecQuantity"
)

setMethod(
    f = "initialize",
    signature = "QuantilePG",
    definition = function(.Object, freqRep, frequencies, levels) {

      .Object@freqRep <- freqRep
      .Object@frequencies <- frequencies
      .Object@levels <- levels

      return(.Object)
    }
)

################################################################################
#' Get values from a quantile periodogram.
#'
#' For vectors \code{frequencies}, \code{levels.1} and \code{levels.2} the
#' values from an \code{object} of type \code{QuantilePG} are returned.
#'
#' Fetching of the periodogram values basically happens by passing
#' \code{frequencies} and the union of \code{levels.1} and \code{levels.2} to
#' \code{\link{getValues}}. Therefore,
#' the parameters \code{frequencies}, \code{levels.1} and \code{levels.1} are
#' expected to be vectors of reals; an error is thrown otherwise. If any of the
#' \code{frequencies}, \code{levels.1} and \code{levels.2} requested is not
#' available from \code{object} a warning is issued. Note that the frequencies
#' are transformed to \eqn{[0,\pi]}{[0,pi]} using \code{\link{frequenciesValidator}}
#' when checking if they are available in \code{object}.
#'
#' The returned array of \code{values} is of dimension \code{[J,K1,K2,B+1]},
#' where \code{J=length(frequencies)}, \code{K1=length(levels.1)},
#' \code{K2=length(levels.2))}, and \code{B} denotes the
#' value stored in slot \code{B} of \code{freqRep} that's a slot of \code{object}.
#' At position \code{(j,k1,k2,b)}
#' the returned value is the one corresponding to \code{frequencies[j]},
#' \code{levels.1[k1]} and \code{levels.2[k2]} that are closest to the
#' \code{frequencies}, \code{levels.1} and \code{levels.2}
#' available in \code{object}; \code{\link{closest.pos}} is used to determine
#' what closest to means.
#'
#' @name getValues-QuantilePG
#' @aliases getValues,QuantilePG-method
#'
#' @keywords Access-functions
#'
#' @param object \code{QuantilePG} of which to get the values
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
#' @examples
#' Y        <- rnorm(32)
#' freq     <- 2*pi*c(0:31)/32
#' levels   <- c(0.25,0.5,0.75)
#' qPG      <- quantilePG(Y, levels.1=levels)
#' V.all    <- getValues(qPG)
#' V.coarse <- getValues(qPG, frequencies = 2*pi*c(0:15)/16)
#' V.fine   <- getValues(qPG, frequencies = 2*pi*c(0:63)/64)
#' V.part   <- getValues(qPG, frequencies = 2*pi*c(0:16)/32,
#'                            levels.1 = c(0.25), levels.2 = c(0.5,0.75))
################################################################################
setMethod(
    f = "getValues",
    signature = signature(
        object="QuantilePG"),
    definition = function(object,
        frequencies=2*pi*(0:(lenTS(object@freqRep@Y)-1))/lenTS(object@freqRep@Y),
        levels.1=getLevels(object,1),
        levels.2=getLevels(object,2),
        d1 = 1:(dim(object@freqRep@Y)[2]),
        d2 = 1:(dim(object@freqRep@Y)[2])) {

    # workaround: default values don't seem to work for generic functions?
    if (!hasArg(frequencies)) {
      frequencies <- 2*pi*(0:(lenTS(object@freqRep@Y)-1))/lenTS(object@freqRep@Y)
    }
    if (!hasArg(levels.1)) {
      levels.1 <- object@levels[[1]]
    }
    if (!hasArg(levels.2)) {
      levels.2 <- object@levels[[2]]
    }
    if (!hasArg(d1)) {
      d1 <- 1:(dim(object@freqRep@Y)[2])
    }
    if (!hasArg(d2)) {
      d2 <- 1:(dim(object@freqRep@Y)[2])
    }
    # end: workaround

    fR <- object@freqRep

    N <- dim(fR@Y)[1]
    D <- dim(fR@Y)[2]
    D1 <- length(d1)
    D2 <- length(d2)
      
    J <- length(frequencies)
    B <- fR@B
    K1 <- length(levels.1)
    K2 <- length(levels.2)

    levels.all <- union(levels.1, levels.2)
    K <- length(levels.all)

    P1 <- match(levels.1,levels.all)
    P2 <- match(levels.2,levels.all)
    x <- getValues(fR, frequencies, levels.all)

    if (D == 1) {
      A <- apply(x, c(1,3), function(x){outer(x[P1],Conj(x[P2]),"*")})
      values <- aperm(array(A, dim=c(K1,K2,J,B+1)), perm=c(3,1,2,4))  
    } else {
      A <- apply(x, c(1,4), function(x){outer(x[,P1],Conj(x[,P2]),"*")})
      values <- aperm(array(A, dim=c(D,K1,D,K2,J,B+1)), perm=c(5,1,2,3,4,6))  
    }
    
    if (D1 == 1 && D2 == 1) {
      final.dim.res <- c(J, K1, K2, object@freqRep@B+1)
    } else {
      final.dim.res <- c(J, D1, K1, D2, K2, object@freqRep@B+1)
    }
    
    values <- array(values, dim=final.dim.res)
    

    return(1/(2*pi*N) * values)
  }
)

################################################################################
#' Get associated \code{\link{FreqRep}} from a
#' \code{\link{QuantilePG}}.
#'
#' @name getFreqRep-QuantilePG
#' @aliases getFreqRep,QuantilePG-method
#'
#' @keywords Access-association-functions
#'
#' @param object \code{QuantilePG} from which to get the
#'                \code{\link{FreqRep}}.
#' @return Returns the \code{\link{FreqRep}} object associated.
################################################################################
setMethod(f = "getFreqRep",
    signature = "QuantilePG",
    definition = function(object) {
      return(object@freqRep)
    }
)

################################################################################
#' Create an instance of the \code{\link{QuantilePG}} class.
#'
#' The parameter \code{type.boot} can be set to choose a block bootstrapping
#' procedure. If \code{"none"} is chosen, a moving blocks bootstrap with
#' \code{l=length(Y)} and  \code{N=length(Y)} would be done. Note that in that
#' case one would also chose \code{B=0} which means that \code{getPositions}
#' would never be called. If \code{B>0} then each bootstrap replication would
#' be the undisturbed time series.
#'
#' @name QuantilePG-constructor
#' @aliases quantilePG
#' @export
#'
#' @keywords Constructors
#'
#' @param Y A \code{vector} of real numbers containing the time series from
#'           which to determine the quantile periodogram or a \code{ts} object
#'           or a \code{zoo} object.
#' @param isRankBased If true the time series is first transformed to pseudo
#'                    data [cf. \code{\link{FreqRep}}].
#' @param levels.1 A vector of length \code{K1} containing the levels \code{x1}
#'                  at which the QuantilePG is to be determined.
#' @param levels.2 A vector of length \code{K2} containing the levels \code{x2}.
#' @param frequencies A vector containing frequencies at which to determine the
#'                     quantile periodogram.
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
#'
#' @return Returns an instance of \code{QuantilePG}.
################################################################################
quantilePG <- function( Y,
                        frequencies=2*pi/lenTS(Y) * 0:(lenTS(Y)-1),
                        levels.1 = 0.5,
                        levels.2=levels.1,
                        isRankBased=TRUE,
                        type=c("clipped","qr"),
                        type.boot = c("none","mbb"),
                        B = 0,
                        l = 0,
                        method = c("br", "fn", "pfn", "fnc", "lasso", "scad"),
                        parallel=FALSE) {

  # Verify if all parameters are valid
  Y <- timeSeriesValidator(Y)

  if (!(is.vector(levels.1)  && is.numeric(levels.1))) {
    stop("'levels.1' needs to be specified as a vector of real numbers")
  }

  if (!(is.vector(levels.2)  && is.numeric(levels.2))) {
    stop("'levels.2' needs to be specified as a vector of real numbers")
  }

  levels.all <- union(levels.1, levels.2)
  K1 <- length(levels.1)
  K2 <- length(levels.2)

  J <- length(frequencies)

  type <- match.arg(type, c("clipped","qr"))[1]
  switch(type,
    "clipped" = {
      freqRep <- clippedFT(Y, frequencies, levels.all, isRankBased, B, l, type.boot)},
    "qr" = {
      freqRep <- qRegEstimator(Y, frequencies, levels.all, isRankBased, B, l, type.boot, method, parallel)}
  )

  obj <- new(
      Class = "QuantilePG",
      freqRep = freqRep,
      frequencies = frequencies,
      levels = list(levels.1, levels.2)
  )

  return(obj)
}

################################################################################
#' Plot the values of the \code{\link{QuantilePG}}.
#'
#' Creates a \code{K} x \code{K} plot depicting a quantile periodogram.
#' Optionally, a simulated copula spectral density can
#' be displayed.
#' In each of the subplots either the real part (on and below the diagonal;
#' i. e., \eqn{\tau_1 \leq \tau_2}{tau1 <= tau2}) or the imaginary parts
#' (above the diagonal; i. e., \eqn{\tau_1 > \tau_2}{tau1 > tau2}) of
#' \itemize{
#'   \item the quantile periodogram (black line),
#'   \item a simulated quantile spectral density (red line),
#' }
#' for the combination of levels \eqn{\tau_1}{tau1} and \eqn{\tau_2}{tau2}
#' denoted on the left and bottom margin of the plot are displayed.
#' 
#' Currently, only the plot for the first component is shown.
#'
#' @name plot-QuantilePG
#' @aliases plot,QuantilePG,ANY-method
#' @export
#'
#' @importFrom abind abind
#'
#' @param x  The \code{\link{QuantilePG}} object to plot
#' @param qsd  a \code{\link{QuantileSD}} object; will be plotted if not
#'              missing.
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
#'                    plotted; default is all available frequencies but 0; if 0 is the
#'                     only available frequency, then only 0 will be used.
#' @param levels a set of levels for which the values are to be plotted.
#'
#' @return Returns the plot described in the Description section.
################################################################################

setMethod(f = "plot",
    signature = "QuantilePG",
    definition = function(x, qsd,
        ratio = 3/2, widthlab = lcm(1), xlab = expression(omega/2*pi), ylab = NULL,
        type.scaling = c("individual", "real-imaginary", "all"),
        frequencies=x@frequencies[-which(x@frequencies == 0)],
        levels=intersect(x@levels[[1]], x@levels[[2]])) {

    def.par <- par(no.readonly = TRUE) # save default, for resetting...

    # workaround: default values don't seem to work for generic functions?
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
    if (!hasArg(frequencies)) {
      if (length(x@frequencies[-which(x@frequencies == 0)]) > 0) {
        frequencies <- x@frequencies[-which(x@frequencies == 0)]
      } else {
        frequencies <- x@frequencies
      }
    }
    if (!hasArg(levels)) {
      levels <- intersect(x@levels[[1]], x@levels[[2]])
    }
    # end: workaround

    if (length(levels) == 0) {
      stop("There has to be at least one level to plot.")
    }

tryCatch({
    N <- length(x@freqRep@Y)
    K <- length(levels)
    # TODO: make a choice of dimensions possible.
    values <- getValues(x, frequencies = frequencies,
        levels.1=levels, levels.2=levels, d1=1, d2=1)

    if (hasArg(qsd)) {
      j.min <- round(min(frequencies*2^8/(2*pi)))
      j.max <- round(max(frequencies*2^8/(2*pi)))
      freq.csd <- 2*pi*(j.min:j.max)/2^8
      # TODO: make a choice of dimensions possible.
      csd <- getValues(qsd, frequencies = freq.csd,
          levels.1=levels, levels.2=levels, d1=1, d2=1)
    }

    X <- frequencies/(2*pi)

    allVals <- array(values[,,,1], dim=c(length(X), K, K))

    if (hasArg(qsd)) {
      allVals <- abind(allVals, csd, along=1)
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
#        allVals <- values[,i1,i2,1]
#        if (hasArg(sCSD)) {
#          allVals <- c(allVals,csd[,i1,i2])
#        }
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
          plot(x=frequencies/(2*pi), y=Re(values[,i1,i2,1]),
              xlim=c(min(X), max(X)), ylim=c(y.min, y.max),
              #ylim=c(min(Re(allVals)), max(Re(allVals))),
              type="l", xlab="", ylab="")
          if (hasArg(qsd)) {
            lines(x=freq.csd/(2*pi), y=Re(csd[,i1,i2]), col="red")
          }
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
          plot(x=frequencies/(2*pi), y=Im(values[,i1,i2,1]),
              xlim=c(min(X), max(X)), ylim=c(y.min, y.max),
              #ylim=c(min(Im(allVals)), max(Im(allVals))),
              type="l", xlab="", ylab="")
          if (hasArg(qsd)) {
            lines(x=freq.csd/(2*pi), y=Im(csd[,i1,i2]), col="red")
          }
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
