#' @include precintcon.ci.per.year.analysis.r
NULL

#' @name ci.per.year
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com}
#' @title Concentration Index per Year
#' @description Calculates the Concentration Index (CI) per year on a daily 
#' precipitation serie.
#' @aliases precintcon.ci.per.year.analysis ci.per.year
#' @param object a daily precipitation serie.
#' @param interval the interval in millimeters applied for calculating the 
#' concentration index. (Default value: 1) 
#' @usage ci.per.year(object, interval = 1)
#' @return A data.frame (precintcon.ci) containing the following variables:
#' \itemize{
#' \item \code{year} is the year.
#' \item \code{ci} is the concentration index.
#' }
#' @seealso 
#' \code{\link{pplot.lorenz}}
#' \code{\link{read.data}}
#' \code{\link{as.daily}}
#' @examples
#' ##
#' # Loading the daily precipitation serie.
#' data(daily)
#' 
#' ##
#' # Performing the Concentration Index Analysis
#' ci.per.year(daily, interval = 1)
#' @keywords concentration index precipitation
#' @export
ci.per.year <- function(object, interval = 1.0) {
   return(precintcon.ci.per.year.analysis(object, interval))  
}