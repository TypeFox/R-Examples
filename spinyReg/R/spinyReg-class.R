##' Class "spinyreg"
##'
##' Class of object returned by the \pkg{spinyreg} function.
##'
##' @section Slots: \describe{
##'
##' \item{\code{coefficients}:}{numeric vector of coefficients with
##' respect to the original input. Contains the intercept if the model
##' owns any.}
##'
##' \item{\code{alpha}:}{numeric scalar.}
##'
##' \item{\code{gamma}:}{numeric scalar.}
##'
##' \item{\code{normx}:}{Vector (class \code{"numeric"}) containing the
##' square root of the sum of squares of each column of the design
##' matrix.}
##'
##' \item{\code{residuals}:}{Vector of residuals.}
##'
##' \item{\code{r.squared}:}{scalar giving the
##' coefficient of determination.}
##'
##' \item{\code{fitted}:}{Vector of fitted values.}
##'
##' \item{\code{monitoring}:}{List (class \code{"list"}) which
##' contains various indicators dealing with the optimization
##' process.}
##'
##' \item{\code{intercept}:}{Logical which indicates if a intercept is
##' included in the model.}}
##'
##' @section Methods:
##' This class comes with the usual \code{predict(object, newx, ...)},
##' \code{fitted(object, ...)}, \code{residuals(object, ...)}, \code{coefficients(object, ...)},
##' \code{print(object, ...)} and \code{show(object)} generic (undocumented) methods.
##'
##' @aliases fitted,spinyreg-method predict,spinyreg-method
##' print,spinyreg-method show,spinyreg-method residuals,spinyreg-method
##' coefficients,spinyreg-method
##'
##' @docType class
##'
##' @keywords class
##'
##' @name spinyreg-class
##' @rdname spinyreg-class
##'
##' @exportClass spinyreg
##' @exportMethod fitted
##' @exportMethod coefficients
##' @exportMethod residuals
##' @exportMethod predict
##' @exportMethod print
##' @exportMethod show
##'
##' @importFrom stats fitted predict residuals
##' @import methods
setClass("spinyreg",
  representation = representation(
     coefficients = "numeric",
     alpha        = "numeric",
     gamma        = "numeric",
     normx        = "numeric",
     fitted       = "numeric",
     residuals    = "numeric",
     r.squared    = "numeric",
     intercept    = "logical",
     monitoring   = "list")
)

setMethod("print", "spinyreg", definition =
   function(x, ...) {
     ncoef <- length(x@coefficients)
     cat("Linear regression with spinyReg.\n")
     if (x@intercept) {
       cat("- number of coefficients:", ncoef - 1,"+ intercept\n")
     } else {
       cat("- number of coefficients:", ncoef,"(no intercept)\n")
     }

     cat("- alpha parameter : ", format(x@alpha, digits = 3),"\n", sep="")
     cat("- gamma parameter : ", format(x@gamma, digits = 3),"\n", sep="")

     invisible(x)
   }
)

setMethod("show", "spinyreg", definition =
   function(object) {print(object)}
)

setMethod("coefficients", "spinyreg", definition =
   function(object, ...) {
     return(object@coefficients)
   }
)

setMethod("fitted", "spinyreg", definition =
   function(object, ...) {
     return(object@fitted)
   }
)

setMethod("predict", "spinyreg", definition =
   function (object, newx=NULL, ...)  {
     if (is.null(newx)) {
       return(object@fitted)
     } else {
       if (object@intercept) {
         cbind(1,newx) %*% object@coefficients
       } else {
         newx %*% object@coefficients
       }
     }
   }
)

setMethod("residuals", "spinyreg", definition =
   function(object, newx=NULL, newy=NULL, ...) {
     if (is.null(newx) | is.null(newy)) {
       return(object@residuals)
     } else {
       return(newy - predict(object, newx))
     }
   }
)
