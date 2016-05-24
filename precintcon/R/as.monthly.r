#' @include as.precintcon.monthly.r
NULL

#' @name as.monthly
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com}
#' @aliases as.precintcon.monthly as.monthly
#' @title Convert a daily precipitation serie to a monthly serie
#' @description Converts a daily precipitation serie to a monthly serie.
#' @usage as.monthly(object)
#' @param object a precintcon.daily object or a data.frame containing 
#' 33 or 3 columns
#' @return A data.frame (precintcon.monthly) containing the following variables:
#' \itemize{
#' \item \code{year} is the year.
#' \item \code{month} is the month.
#' \item \code{precipitation} is the precipitation amount in millimeters.
#' }
#' @seealso 
#' \code{\link{pplot.lorenz}}
#' \code{\link{read.data}}
#' @examples 
#' ## Loading the daily precipitation serie.
#' #
#' data(daily)
#' 
#' ## Converting precipitation
#' #
#' as.monthly(daily)
#' @keywords monthly precipitation
#' @export
as.monthly <- function(object) {
   return(as.precintcon.monthly(object))   
}