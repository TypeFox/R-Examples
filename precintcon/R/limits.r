#' @include precintcon.limits.analysis.r
NULL

#' @name limits
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com}
#' @aliases precintcon.limits.analysis limits
#' @title Limits
#' @description Calculates the upper and lower limit on a set of daily or monthly 
#' precipitation series.
#' @usage limits(\dots)
#' @param \dots a set of daily or monthly precipitation series.
#' @return A data.frame containing the following variables:
#' \itemize{
#' \item \code{dataset} is the precipitation serie name.
#' \item \code{max} is the maximum value in the precipitation serie.
#' \item \code{max.date} is the first date of the maximum precipitation serie.
#' \item \code{min} is the minimum value in the precipitation serie.
#' \item \code{min.date} is the first date of the minimum precipitation serie.
#' }
#' @seealso
#' \code{\link{read.data}}
#' \code{\link{as.daily}}
#' @examples 
#' ##
#' # Loading the monthly precipitation serie.
#' data(monthly)
#' 
#' ## 
#' # Performing the Concentration Index Analysis
#' limits(monthly)
#' @keywords limits precipitation
#' @export
limits <- function(...) {
   return(precintcon.limits.analysis(...))
}