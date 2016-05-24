
.subset <- function (x, ...){
  validObject (x)
  x@data <- subset (x@data, ...)
  validObject (x)
  
	x
}


##' subset for hyperSpec object
##'
##' @title subset
##' @name subset
##' @param x hyperSpec object
##' @param ... handed to \code{\link[base]{subset}} (data.frame method)
##' @docType methods
##' @aliases subset subset,hyperSpec-method
##' @return hyperSpec object containing the respective subset of spectra.
##' @author Claudia Beleites
##' @seealso \code{\link[base]{subset}}
##' @export 
setMethod ("subset", signature = signature (x = "hyperSpec"), .subset)
