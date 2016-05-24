#' @include precintcon.rai.analysis.r
NULL

#' @name rai 
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com} 
#' @aliases precintcon.rai.analysis rai 
#' @title Rainfall Anomaly Index 
#' @description Calculates the Rainfall Anomaly Index (RAI) on a daily or 
#' monthly precipitation serie. 
#' @usage rai(object, granularity = "m") 
#' @param object a daily or monthly precipitation serie.
#' @param granularity the granularity applied for calculating the rainfall anomaly index, 
#'   which should be either "m" for monthly granularity or "a" for annual granularity. (Default value: "m")
#' @return A data.frame (precintcon.rai) containing the following variables:
#' \itemize{
#' \item \code{year} is the year.
#' \item \code{month} is the month. This attribute exists only when granularity = "m".
#' \item \code{rai} is the rainfall anomaly index.
#' }
#' @seealso 
#' \code{\link{pplot.rai}}
#' \code{\link{read.data}}
#' @examples 
#' ##
#' # Loading the monthly precipitation serie.
#' data(monthly)
#' 
#' ##
#' # Performing the Rainfall Anomaly Index analysis
#' rai(monthly, granularity = "m")
#' @references Van Rooy, M. P. "A rainfall anomaly index independent of time and space." Notos 14.43 (1965): 6.
#' @keywords rainfall anomaly index precipitation
#' @export 
rai <- function(object, granularity = "m") precintcon.rai.analysis(object, granularity)