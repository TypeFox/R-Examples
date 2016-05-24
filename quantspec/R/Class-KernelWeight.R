#' @include generics.R
#' @include Class-Weight.R
NULL

################################################################################
#' Class for Brillinger-type Kernel weights.
#'
#' \code{KernelWeight} is an S4 class that implements a weighting function by
#' specification of a kernel function \code{W} and a scale parameter \code{bw}.
#'
#' It extends the class \code{\link{Weight}} and writes
#' \deqn{W_N(2\pi (k-1)/N) := \sum_{j \in Z} bw^{-1} W(2\pi bw^{-1} [(k-1)/N + j])}
#' to \code{values[k]} [nested inside \code{env}] for \code{k=1,\ldots,N}.
#' The number \code{length(values)} of Fourier frequencies for which
#' \eqn{W_N} will be evaluated may be set on construction or updated when
#' evoking the method \code{getValues}.
#' To standardize the weights used in the convolution to unity
#' \deqn{W_N^j := \sum_{j \neq s = 0}^{N-1} W_n(2\pi s / N)}
#' is stored to \code{Wnj[s]} for \code{s=1,\ldots,N}, for later usage.
#'
#' @name   KernelWeight-class
#' @aliases KernelWeight
#' @exportClass KernelWeight
#'
#' @keywords S4-classes
#'
#' @slot W a kernel function
#' @slot bw bandwidth
#' @slot env An environment to allow for slots which need to be
#'       accessable in a call-by-reference manner:
#'       \describe{
#'         \item{\code{values}}{A vector storing the weights; see the Details
#'               section.}
#'         \item{\code{Wnj}}{A vector storing the terms used for normalization;
#'               see the Details section.}
#'       }
#'
#' @seealso
#' Examples for implementations of kernels \code{W} can be found at:
#' \code{\link{kernels}}.
#'
#' @references
#' Brillinger, D. R. (1975). \emph{Time Series: Data Analysis and Theory}.
#' Holt, Rinehart and Winston, Inc., New York. [cf. p. 146 f.]
################################################################################
setClass(
    Class = "KernelWeight",
    representation=representation(
        W   = "function",
        bw  = "numeric",
        env = "environment"   # Enviroment 'hosting':
                              #   o values = "numeric", and
                              #   o Wnj    = "numeric".
    ),
    contains = "Weight"
)

setMethod(
    f = "initialize",
    signature = "KernelWeight",
    definition = function(.Object, W, bw, N, descr) {

      .Object@W <- W
      .Object@bw <- bw
      .Object@descr <- descr

      .Object@env <- new.env(parent=emptyenv())
      .Object@env$values <- W(0)/bw
      .Object@env$values <- getValues(.Object, N=N)

      return(.Object)
    }
)

################################################################################
#' Get values from a weight object determined by a kernel function \code{W} and
#' a bandwidth \code{b}.
#'
#' For an \code{object} of type \code{KernelWeight} and an optional integer
#' \code{N} the weights \eqn{W_n}{Wn} are returned as a vector that has
#' \eqn{W_n(2\pi(k-1)/n)}{Wn(2pi(k-1)/n)} at position \code{k}.
#'
#' @name getValues-KernelWeight
#' @aliases getValues,KernelWeight-method
#'
#' @keywords Access-functions
#'
#' @param object \code{KernelWeight} of which to get the values
#' @param N a \code{numeric} specifying the number of equaly spaced Fourier
#'        frequencies from \eqn{[0,2\pi)}{[0,2pi)} for which the weight will
#'        be computed; by default the number \code{N} specified on construction.
#' @return Returns a vector of size \code{N} as described in the Details section.
################################################################################
setMethod(f = "getValues",
    signature = "KernelWeight",
    definition = function(object, N = length(object@env$values)) {
      if (length(object@env$values) != N) {

        bw <- object@bw
        W <- object@W
        Wn <- function(x) {
          start <- ceiling(-bw/2-x/(2*pi))
          end <- floor(bw/2-x/(2*pi))
          result <- 0
          for (j in start:end) {
            result <- result + bw^(-1) * W(bw^(-1)*(x+2*pi*j))
          }
          result
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
#' Get attribute \code{W} (kernel used for smoothing) from a \code{KernelWeight}.
#'
#' @name getW-KernelWeight
#' @aliases getW,KernelWeight-method
#'
#' @keywords Access-functions
#'
#' @param object \code{KernelWeight} from which to get the kernel \code{W}.
#' @return Returns the \code{W} attribute.
################################################################################
setMethod(f = "getW",
    signature = "KernelWeight",
    definition = function(object) {
      return(object@W)
    }
)

################################################################################
#' Get attribute \code{bw} (bandwidth / scaling parameter used for smoothing)
#' from a \code{KernelWeight}.
#'
#' @name getBw-KernelWeight
#' @aliases getBw,KernelWeight-method
#'
#' @keywords Access-functions
#'
#' @param object \code{KernelWeight} from which to get the bandwidth \code{bw}.
#' @return Returns the \code{bw} attribute.
################################################################################
setMethod(f = "getBw",
    signature = "KernelWeight",
    definition = function(object) {
      return(object@bw)
    }
)

################################################################################
#' Get attribute \code{Wnj} from a \code{QSpecQuantity}.
#'
#' If the optional parameter \code{j} is supplied, then only the \code{j}th
#' element(s) of the vector will be returned, the entire vector otherwise.
#'
#' @name getWnj-KernelWeight
#' @aliases getWnj,KernelWeight-method
#'
#' @keywords Access-functions
#'
#' @param object \code{KernelWeight} from which to get the \code{Wnj}.
#' @param j an integer or vector of indices specifying which \code{Wnj[j]}
#'           to return.
#'
#' @return Returns levels attribute, as a vector of real numbers.
#'
#' @examples
#' wgt <- kernelWeight(W=W1, N=2^3, bw=0.7)
#' getWnj(wgt)
#' getWnj(wgt, 2)
#' getWnj(wgt, c(2,7))
################################################################################
setMethod(f = "getWnj",
    signature = "KernelWeight",
    definition = function(object,j) {
      if (missing("j")) {
        return(object@env$Wnj)
      } else {
        is.wholenumber <-
            function(x, tol = .Machine$double.eps^0.5) {abs(x - round(x)) < tol}

        if (!all(is.wholenumber(j) & (j > 0))) {
          stop("j needs to be a vector of positive integers")
        } else {
          return(object@env$Wnj[j])
        }
      }
    }
)


################################################################################
#' Create an instance of the \code{\link{KernelWeight}} class.
#'
#' @name KernelWeight-constructor
#' @aliases kernelWeight
#' @export
#'
#' @keywords Constructors
#'
#' @param W A kernel function
#' @param N Fourier basis; number of grid points in \eqn{[0,2\pi)}{[0,2pi)} for
#'           which the weights will be computed.
#' @param bw bandwidth; if a vector, then a list of weights is returned
#' @param descr a description to be used in some plots
#'
#' @return Returns an instance of \code{KernelWeight}.
#'
#' @seealso \code{\link{kernels}}
#'
#' @examples
#' wgt1 <- kernelWeight(W=W0, N=16, bw=c(0.1,0.3,0.7))
#' print(wgt1)
#' wgt2 <- kernelWeight(W=W1, N=2^8, bw=0.1)
#' plot(wgt2, main="Weights determined from Epanechnikov kernel")
################################################################################
kernelWeight <- function(W = W0,
    N=1, bw = 0.1*N^(-1/5), descr = paste("bw=",round(bw,3),", N=",N,sep="")) {

  res <- c()

  for (i in 1:length(bw)) {
    obj <- new(
        Class = "KernelWeight",
        W, bw[i], N, descr[i]
    )
    res <- c(res,obj)
  }
  if (length(res) == 1) {res <- res[[1]]}
  return(res)
}


setMethod(f="show",
    signature = signature(object="KernelWeight"),
    definition = function(object) {
      cat("\n",class(object),": ",object@descr)
      cat("\n")
      print(getValues(object))
    }
)

################################################################################
#' Plot the values of the \code{\link{KernelWeight}}.
#'
#' Creates a plot visualizing the weights \eqn{W_n(\omega)}
#' [cf. \code{\link{KernelWeight-class}}] that are used to estimate the
#' quantile spectral density.
#'
#' In the plot the values at the frequencies
#' \eqn{2\pi j / N}, \eqn{j=L+1-N,\ldots,L}, \eqn{L:=\lfloor N/2 \rfloor} are
#' shown, where \eqn{N} is the parameter specified on construction of the object
#' or \eqn{N := 3}, if that parameter was smaller than three. A warning is given
#' in the later case.
#'
#' @name plot-KernelWeight
#' @aliases plot,KernelWeight,missing-method
#' @export
#'
#' @param x  The \code{\link{KernelWeight}} to plot.
#' @param y missing arg from the generic; will be ignored.
#' @param xlab label for the x-axis; optional
#' @param ylab label for the y-axis; optional
#' @param main titel (on top) of the plot; optional
#' @param ... optional parameters used for plotting
#'
#' @return Plots the \code{\link{KernelWeight}}.
#'
#' @examples
#' plot(kernelWeight(W1, bw=0.3),
#'        ylab=expression(W[n](x)),
#'        xlab=expression(x),
#'        main="Weights to an Epanechnikov kernel", sub="bw=0.3")
################################################################################
setMethod(f="plot",
    signature = signature(x="KernelWeight", y="missing"),
    definition = function(x, y,
        ylab=expression(W[n](omega)),
        xlab=expression(omega),
        main=x@descr, ...) {
      N <- N.old <- length(getValues(x))
      if (N < 3) {
        warning("N < 3 is not large enough for a nice plot; showing N=128")
        N <- 128
        V <- getValues(x, N=128)
      }
      L <- floor(N/2)
      X <- 2*pi*((L+1-N):L)/N
      Y <- x@env$values[c((L+2):N,1:(L+1))]
      plot(x=X, y=Y, type="l",
          ylab=ylab, xlab=xlab, main=main, ...)
      V <- getValues(x, N=N.old)
    }
)
