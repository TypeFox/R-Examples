#' @include as.precintcon.seasonal.r
NULL

#' @name as.seasonal
#' @aliases as.precintcon.seasonal as.seasonal 
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com} 
#' 
#' @title Converting to seasonal precipitation serie
#' @description Converts a daily or monthly precipitation serie 
#' to a seasonal serie according to meteorological seasons.
#' @details If the serie has no a month of a season, it is waived in convertion, e.g., 
#' if a serie has January and February of 1975, but no December of 1974, the first two months are
#' removed of the resulting serie because the season that depends all of them is not complete.
#' @usage as.seasonal(object) 
#' @param object a precintcon.daily, or precintcon.monthly object or 
#' a data.frame containing 33 or 3 columns.
#' @return A data.frame (precintcon.seasonal) containing the following variables:
#' \itemize{
#' \item \code{year} is the year.
#' \item \code{season} is the season.
#' \item \code{precipitation} is the precipitation amount in millimeters.
#' }
#' @seealso 
#' \code{\link{pplot.lorenz}}
#' \code{\link{read.data}}
#' @examples 
#' ##
#' # Loading the daily precipitation serie.
#' data(daily)
#' 
#' ##
#' # Converting precipitation
#' as.seasonal(daily)
#' @keywords seasonal precipitation
#' @export 
as.seasonal <- function(object) as.precintcon.seasonal(object)