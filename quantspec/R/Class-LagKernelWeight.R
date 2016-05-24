#' @include generics.R
#' @include Class-Weight.R
#' @include kernels.R
NULL

################################################################################
#' Class for lag window generators
#'
#' \code{LagKernelWeight} is an S4 class that implements a weighting function by
#' specification of a kernel function \code{W} and a scale parameter \code{bw}.
#'
#' It extends the class \code{\link{Weight}} and writes
#' \deqn{W_N(x[k]) := W(x[k]/bw)}
#' to \code{values[k]} [nested inside \code{env}] for \code{k=1,\ldots,length(x)}.
#' The points \code{x} where \eqn{W} is evaluated may be set on construction or updated when
#' evoking the method \code{getValues}.
#'
#' @name   LagKernelWeight-class
#' @aliases LagKernelWeight
#' @exportClass LagKernelWeight
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
#'       }
#'
#' @seealso
#' Examples for implementations of kernels \code{W} can be found at:
#' \code{\link{kernels}}.
#'
################################################################################
setClass(
  Class = "LagKernelWeight",
  representation=representation(
    W   = "function",
    bw  = "numeric",
    env = "environment"   # Enviroment 'hosting':
    #   o values = "numeric"
  ),
  contains = "Weight"
)

setMethod(
  f = "initialize",
  signature = "LagKernelWeight",
  definition = function(.Object, W, bw, K,  descr) {
    
    .Object@W <- W
    .Object@bw <- bw
    .Object@descr <- descr
    
    .Object@env <- new.env(parent=emptyenv())
    .Object@env$values <- W(0)/bw
    .Object@env$values <- getValues(.Object, K = K)
    
    return(.Object)
  }
)

################################################################################
#' Get values from a weight object determined by a kernel function \code{W} and
#' a bandwidth \code{bw}.
#'
#' For an \code{object} of type \code{LagKernelWeight} and an optional integer
#' \code{K} the weights \eqn{W_k}{Wk} are returned as a vector that has
#' \eqn{W_k((k-1)/bw)}{Wk((k-1)/bw)} at position \code{k}.
#'
#' @name getValues-LagKernelWeight
#' @aliases getValues,LagKernelWeight-method
#'
#' @keywords Access-functions
#'
#' @param object \code{LagKernelWeight} of which to get the values
#' @param K a \code{numeric} that determines the largest lag. The weight will be
#'        computed for the \eqn{K} integers \eqn{0:(K-1)}; by default the number
#'        \code{K} specified on construction.
#' @return Returns a vector of size \code{K} as described in the Details section.
################################################################################
setMethod(f = "getValues",
          signature = "LagKernelWeight",
          definition = function(object, K = length(object@env$values)) {
            if (length(object@env$values) != K) {
            	Wn = object@W
		bw = object@bw
		object@env$values <- Wn((0:(K-1))/bw)
            }
            return(object@env$values)
          }
)
################################################################################
#' Get attribute \code{W} (kernel used for smoothing) from a \code{LagKernelWeight}.
#'
#' @name getW-LagKernelWeight
#' @aliases getW,LagKernelWeight-method
#'
#' @keywords Access-functions
#'
#' @param object \code{LagKernelWeight} from which to get the kernel \code{W}.
#' @return Returns the \code{W} attribute.
################################################################################
setMethod(f = "getW",
          signature = "LagKernelWeight",
          definition = function(object) {
            return(object@W)
          }
)

################################################################################
#' Get attribute \code{bw} (bandwidth / scaling parameter used for smoothing)
#' from a \code{LagKernelWeight}.
#'
#' @name getBw-LagKernelWeight
#' @aliases getBw,LagKernelWeight-method
#'
#' @keywords Access-functions
#'
#' @param object \code{LagKernelWeight} from which to get the bandwidth \code{bw}.
#' @return Returns the \code{bw} attribute.
################################################################################
setMethod(f = "getBw",
          signature = "LagKernelWeight",
          definition = function(object) {
            return(object@bw)
          }
)


################################################################################
#' Create an instance of the \code{\link{LagKernelWeight}} class.
#'
#' @name LagKernelWeight-constructor
#' @aliases lagKernelWeight
#' @export
#'
#' @keywords Constructors
#'
#' @param W A kernel function
#' @param K a \code{numeric} that determines the largest lag. The weight will be
#'        computed for the \eqn{K} integers \eqn{0:(K-1)}; by default the number
#'        \code{K} specified on construction.
#' @param bw bandwidth
#' @param descr a description to be used in some plots
#'
#' @return Returns an instance of \code{LagKernelWeight}.
#'
#' @seealso \code{\link{kernels}}
#'
#' @examples
#' wgt1 <- lagKernelWeight(W=WParzen, K=20, bw=10)
#' print(wgt1)
################################################################################
lagKernelWeight <- function(W = WParzen,  bw = K/2, K = 10, descr = paste("bw=",bw,", K=",K,sep="")) {
  
  #TODO test if parameters are okay
  Kernel = new(Class = "LagKernelWeight", W=W,K=K,bw=bw,descr = descr)
  return(Kernel)
}


setMethod(f="show",
          signature = signature(object="LagKernelWeight"),
          definition = function(object) {
            cat("\n",class(object),": ",object@descr)
            cat("\n")
            print(getValues(object))
          }
)

################################################################################
#' Plot the values of the \code{\link{LagKernelWeight}}.
#'
#' Creates a plot visualizing the weights \eqn{W_n(k}
#' [cf. \code{\link{LagKernelWeight-class}}] that are used to estimate the
#' quantile spectral density.
#'
#' In the plot the values at the points \eqn{k/bw}{k/bw} with 
#' \eqn{k \in \{-K,\dots,K\}}{(-K:k:K)} are shown.
#'
#' @name plot-LagKernelWeight
#' @aliases plot,LagKernelWeight,missing-method
#' @export
#'
#' @param x  The \code{\link{LagKernelWeight}} to plot.
#' @param y missing arg from the generic; will be ignored.
#' @param xlab label for the x-axis; optional
#' @param ylab label for the y-axis; optional
#' @param main titel (on top) of the plot; optional
#' @param ... optional parameters used for plotting
#'
#' @return Plots the \code{\link{LagKernelWeight}}.
#'
#' @examples
#' plot(lagKernelWeight(WParzen, bw=10, K = 20),
#'        ylab=expression(W[n](x)),
#'        xlab=expression(x),
#'        main="Weights to the Parzen Window")
################################################################################
setMethod(f="plot",
          signature = signature(x="LagKernelWeight", y="missing"),
          definition = function(x, y,
                                ylab=expression(W[n](k)),
                                xlab=expression(k),
                                main=x@descr, ...) {
            K <- K.old <- length(getValues(x))
            if (K < 3) {
              warning("K < 3 is not large enough for a nice plot; showing N=128")
              N <- 128
              V <- getValues(x, K=128)
            }
            X = -(K-1):(K-1)
            V = getValues(x)
	    Y = c(V[K:1],V[-1])
            plot(x=X, y=Y, type="l",
                 ylab=ylab, xlab=xlab, main=main, ...)
            V <- getValues(x, K=K.old)
          }
)
