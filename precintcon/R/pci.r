#' @include precintcon.pci.analysis.r
NULL

#' @name pci
#' @aliases precintcon.pci.analysis pci
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com}
#' 
#' @title Precipitation Concentration Index
#' @description Calculates the Precipitation Concentration Index (PCI) on a 
#' daily or monthly precipitation serie.
#' @usage pci(object)
#' @param object a daily or monthly precipitation serie.
#' @return A data.frame containing the following variables:
#' \itemize{
#' \item \code{year} is the year.
#' \item \code{pci} is the precipitation concentration index.
#' }
#' @seealso 
#' \code{\link{pplot.pci}}
#' \code{\link{read.data}}
#' \code{\link{as.daily}}
#' \code{\link{as.monthly}}
#' @examples 
#' ##
#' # Loading the monthly precipitation serie.
#' data(monthly)
#' 
#' ## 
#' # Performing the Precipitation Concentration Index analysis
#' pci(monthly)
#' @keywords precipitation concentration index precipitation
#' @export
pci <- function(object) {
   return(precintcon.pci.analysis(object))  
}