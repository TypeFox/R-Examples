#' @title Generic methods for objects with class \code{'spec'}
#'
#' @details Objects with class \code{'spec'} are simply lists with spectral estimates and parameters 
#' \code{as.data.frame} converts the list into a \code{'data.frame'} with individual
#' columns for the frequency, PSD, and taper vectors; 
#' all other information will be retained as an attribute.
#'
#' @name spec-methods
#' @author A.J. Barbour
#' @rdname spec-methods
#' @docType methods
#' 
#' @param x a \code{'spec'} object
#' @param y optional coordinate vector for the y-axis
#' @param type character; the type of plot
#' @param ... optional arguments
#' 
#' @example inst/Examples/rdex_spec.R
NULL

#' @rdname spec-methods
#' @export
as.list.spec <- function(x, ...){
  class(x) <- 'list'
  return(x)
}

#' @rdname spec-methods
#' @export
as.spec <- function(x, ...) UseMethod("as.spec")

#' @rdname spec-methods
#' @export
as.spec.amt <- function(x, ...){
  class(x) <- 'spec'
  return(x)
}

#' @rdname spec-methods
#' @aliases lines.spec
#' @export
lines.spec <- function(x, y=NULL, type = 'l', ...){
  plot(x, add=TRUE, type=type, ...)
}

#' @rdname spec-methods
#' @aliases as.data.frame.spec
#' @export
as.data.frame.spec <- function(x, ...){
  
  # get the meat-n-potatoes
  xdf <- data.frame(freq=x[['freq']], spec=x[['spec']], taper = x[['taper']])
  
  # keep all else as an attribute -- is there a better way?
  attr(xdf, "coh")       <- x[['coh']]
  attr(xdf, "phase")     <- x[['phase']]
  attr(xdf, "kernel")    <- x[['kernel']]
  attr(xdf, "df")        <- x[['df']]
  attr(xdf, "bandwidth") <- x[['bandwidth']]
  attr(xdf, "n.used")    <- x[['n.used']]
  attr(xdf, "orig.n")    <- x[['orig.n']]
  attr(xdf, "series")    <- x[['series']]
  attr(xdf, "snames")    <- x[['snames']]
  attr(xdf, "method")    <- x[['method']]
  attr(xdf, "pad")       <- x[['pad']]
  attr(xdf, "detrend")   <- x[['detrend']]
  attr(xdf, "demean")    <- x[['demean']]
  
  return(xdf)
}

#' @rdname spec-methods
#' @export
data.frame.spec <- as.data.frame.spec
