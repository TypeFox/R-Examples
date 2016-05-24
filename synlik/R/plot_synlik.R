##########
#### Method to plot the object
##########

.plot.synlik <- function(x, ...)
{
  if(!is(x, "synlik")) stop("object has to be of class \"synlik\" ")
  
  tmpPlot <- x@plotFun
  
  if( !is.null(tmpPlot) ){ 
    tmpPlot(x@data, ...)
  } else {
    e <- tryCatch(plot(x@data, ...), error = function(e) e)
    if("error" %in% class(e)) message("I don't know how to plot object@data, put the right plotting function in object@plotFun.")
  }
  
  return( invisible() )
}

#' @name plot-synlik
#' 
#' @title Method for plotting an object of class \code{synlik}.
#'
#' @description It basically calls the slot \code{object@@plotFun} with input \code{object@@data}, if it has been provided by the user. 
#' Otherwise it tries to use the \code{plot(x = object@@data, y, ...)} generic.
#'
#' @param x An object of class \code{synlik}.
#' @param ... additional arguments to be passed to \code{object@@plotFun}.
#'
#' @return NULL
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>  
#' @seealso \code{\link{synlik-class}}, \code{\link{plot}}.
#' @aliases plot,synlik,missing-method
#' @method plot synlik missing
#' @examples
#' data(ricker_sl)
#' 
#' # Using ricker_sl@@plotFun
#' plot(ricker_sl)
#' 
#' # Using generic plot, doesn't work well because object@@data is a matrix. 
#' ricker_sl@@plotFun <- NULL
#' plot(ricker_sl)
#' 
#' @rdname plot-synlik
setMethod("plot",
          signature = signature(x = "synlik", y = "missing"),
          definition = .plot.synlik)