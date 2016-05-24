#' @include generics.R
#' @include Class-Weight.R
NULL

################################################################################
#' Class for weights to estimate integrated spectral density kernels.
#'
#' \code{SpecDistrWeight} is an S4 class that implements a weighting function given
#' by
#' \deqn{W_n(\alpha) := I\{\alpha \leq 0\}}{Wn(a) := I{a <= 0}}.
#'
#' At position \code{k} the value \eqn{W_n(2\pi (k-1)/n}{Wn(2pi(k-1)/n)} is
#' stored [in a vector \code{values} nested inside \code{env}] for \code{k=1,\ldots,T}.
#' The number \code{length(values)} of Fourier frequencies for which
#' \eqn{W_n}{Wn} will be evaluated may be set on construction or updated when
#' evoking the method \code{getValues}.
#'
#' @name   SpecDistrWeight-class
#' @aliases SpecDistrWeight
#' @exportClass SpecDistrWeight
#'
#' @keywords S4-classes
################################################################################
setClass(
    Class = "SpecDistrWeight",
    contains = "Weight"
)

setMethod(
    f = "initialize",
    signature = "SpecDistrWeight",
    definition = function(.Object, descr) {
      #cat("~~~ SpecDistrWeight: initializator ~~~ \n")

      .Object@descr <- descr

      return(.Object)
    }
)

################################################################################
#' Get values from a weight object of type \code{\link{SpecDistrWeight}}
#'
#' For an \code{object} of type \code{SpecDistrWeight} and an optional integer
#' \code{N} the weights \eqn{W_n}{Wn} are returned as a vector that has
#' \eqn{W_n(2\pi(k-1)/n)}{Wn(2pi(k-1)/n)} at position \code{k}.
#'
#' @name getValues-SpecDistrWeight
#' @aliases getValues,SpecDistrWeight-method
#'
#' @param object \code{SpecDistrWeight} of which to get the values
#' @param N a \code{numeric} specifying the number of equaly spaced Fourier
#'            frequencies from \eqn{[0,2\pi)}{[0,2pi)} for which the weight will
#'            be computed; by default the number \code{N} specified on construction.
#'
#' @return Returns a vector of size \code{N} as described in the Description
#'         section.
################################################################################
setMethod(f = "getValues",
    signature = "SpecDistrWeight",
    definition = function(object, N = length(object@env$values)) {

    if (length(object@env$values) != N) {

      Wn <- function(x) {
        return((x <= 0))
      }

      v <- unlist(lapply(2*pi * (0:(N-1))/N,Wn))
      object@env$values <- v
      object@env$Wnj <- (2*pi/N) * (sum(v) - v[c(2:N,1)])
      object@env$Wnj[object@env$Wnj == 0] <- 1
    }
    return(object@env$values)
  }
)

################################################################################
#' Create an instance of the \code{\link{SpecDistrWeight}} class.
#'
#' @name SpecDistrWeight-constructor
#' @aliases specDistrWeight
#' @export
#'
#' @keywords Constructors
#'
#' @param descr a description for the weight object
#'
#' @return an instance of \code{SpecDistrWeight}.
#'
#' @examples
#' wgt <- specDistrWeight()
################################################################################
specDistrWeight <- function(descr = "Spectral Distribution Weights") {

  obj <- new(
      Class = "SpecDistrWeight",
      descr
  )

  return(obj)
}

setMethod(f="show",
    signature = signature(object="SpecDistrWeight"),
    definition = function(object) {
      cat("\n",class(object),": ",object@descr)
    }
)

################################################################################
#' Plot the values of the \code{\link{SpecDistrWeight}}.
#'
#' Creates a plot visualizing the weights \eqn{W_n(\omega)}
#' [cf. \code{\link{SpecDistrWeight-class}}] that are used to estimate the
#' integrated quantile spectral density.
#'
#' In the plot the values at the frequencies
#' \eqn{2\pi j / 128}, \eqn{j=-63,\ldots,64} are shown.
#'
#' @name plot-SpecDistrWeight
#' @aliases plot,SpecDistrWeight,missing-method
#' @export
#'
#' @param x The \code{\link{SpecDistrWeight}} to plot.
#' @param y missing arg from the generic; will be ignored.
#' @param xlab label for the x-axis; optional
#' @param ylab label for the y-axis; optional
#' @param main titel (on top) of the plot; optional
#' @param ... optional parameters used for plotting
#'
#' @return Plots the \code{\link{SpecDistrWeight}}.
#'
#' @examples
#' plot(specDistrWeight(),
#'        ylab=expression(W[n](x)),
#'        xlab=expression(x))
################################################################################
setMethod(f="plot",
    signature = signature(x="SpecDistrWeight", y="missing"),
    definition = function(x, y,
        ylab=expression(W[n](omega)),
        xlab=expression(omega),
        main=x@descr, ...) {
      N <- 128
      L <- floor(N/2)
      X <- 2*pi*((L+1-N):L)/N
      Y <- (X<=0)
      plot(x=X, y=Y, type="l",
          ylab=ylab, xlab=xlab, main=main, ...)
    }
)
