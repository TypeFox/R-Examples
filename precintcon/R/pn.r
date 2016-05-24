#' @include precintcon.pn.analysis.r
NULL

#' @name pn
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com} 
#' @aliases precintcon.pn.analysis pn 
#' @title Percentage of Normal 
#' @description Calculates the Percentage of Normal (PN) on a daily or 
#' monthly precipitation serie. 
#' @usage pn(object, interval = 30, scale = "a") 
#' @param object a daily or monthly precipitation serie.
#' @param interval the number of months applied for calculating the percentage of 
#' normal.
#' @param scale the scale used for calculating the percentage of normal, 
#' which should be either "w" for weak (not supported yet), 
#' "m" for month, "s" for season, or "d" for decades.
#' @return A data.frame (precintcon.pn) containing the following variables:
#' \itemize{ 
#' \item \code{year} is the year.
#' \item \code{month} is the month. It exists only whether scale = "m".
#' \item \code{pn} is the percentage of normal.
#' }
#' @seealso 
#' \code{\link{pplot.pn}}
#' \code{\link{read.data}}
#' \code{\link{as.daily}}
#' \code{\link{as.monthly}}
#' @examples 
#' ##
#' # Loading the monthly precipitation serie.
#' data(monthly)
#' 
#' ##
#' # Performing the Percentage of Normal analysis
#' pn(monthly)
#' @keywords precipitation percent of normal
#' @export
pn <- function(
   object, 
   interval = 30, 
   scale    = "a"
) {
   return(precintcon.pn.analysis(object, interval, scale))
}