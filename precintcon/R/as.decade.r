#' @include as.precintcon.decade.r
NULL

#' @name as.decade
#' @aliases as.precintcon.decade as.decade
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com}
#' 
#' @title Converting a precipitation serie to a decade serie
#' @description Converts a daily, monthly or annual 
#' precipitation serie to a decade serie.
#' @details It excludes no complete decades for converting the serie, e.g., a serie starting in 1977 and
#' finishing in 2008 will have the year 1977 to 1979 and 2000 to 2008 excluded, resulting into a serie of
#' the years 1980 and 1990.
#' @usage as.decade(object)
#' @param object a precintcon.daily, precintcon.monthly, precintcon.seasonal,
#' precintcon.annual object or a data.frame containing 33 or 3 columns.
#' @return A data.frame (precintcon.decade) containing the following variables:
#' \itemize{
#'	\item \code{year} is the year.
#' \item \code{precipitation} is the decade's precipitation in millimeters.   
#' }
#' @seealso 
#' \code{\link{as.precintcon.annual}} 
#' \code{\link{as.precintcon.seasonal}} 
#' \code{\link{as.precintcon.monthly}}
#' \code{\link{as.precintcon.daily}}
#' @examples 
#' ##
#' # Loading the daily precipitation serie.
#' data(daily)
#' 
#' ##
#' # Converting precipitation
#' as.decade(daily)
#' @keywords decade precipitation
#' @export
as.decade <- function(object) as.precintcon.decade(object)