#'Get PSR data.
#'
#'@description Extract Plant Summary Report from a YplantQMC simulation.
#'
#'Simple wrapper: \code{x$psrdata} is the same as \code{psrdata(x).}
#'
#'@param x Object returned by \code{\link{YplantDay}}
#'@param \dots Further arguments ignored
#'@author Remko Duursma
#'@keywords misc
#'@export
psrdata <- function(x,...){
	if(!inherits(x,"yplantsim"))stop("Need object of class 'yplantsim'.")
	
	x$psrdata

}
