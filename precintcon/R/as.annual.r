#' @include as.precintcon.annual.r
NULL

#' @name as.annual
#' @aliases as.precintcon.annual as.annual
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com}
#' 
#' @title Converting to seasonal precipitation serie.
#' @description Converts a daily, monthly, or seasonal precipitation serie to an 
#' annual serie.
#' @usage as.annual(object)
#' @param object a precintcon.daily, precintcon.monthly, 
#' or precintcon.seasonal object or a data.frame containing 33 or 3 columns.
#' @return A data.frame (precintcon.annual) containing the following variables:
#' \itemize{
#' \item \code{year} is the year.
#' \item \code{precipitation} is the precipitation amount in millimeters.
#' }
#' @seealso 
#' \code{\link{as.daily}}
#' \code{\link{as.monthly}}
#' \code{\link{as.seasonal}}
#' \code{\link{pplot.lorenz}}
#' \code{\link{read.data}}
#' @keywords annual precipitation
#' @examples
#' ##
#' # Loading the daily precipitation serie.
#' data(daily)
#' 
#' ##
#' # Converting precipitation
#' as.annual(daily)
#' @export
as.annual <- function(object) {
   return(as.precintcon.annual(object))
}