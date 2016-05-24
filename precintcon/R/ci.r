#' @include precintcon.ci.analysis.r
NULL

#' @name ci
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com}
#' @aliases precintcon.ci.analysis ci
#' @title Concentration Index
#' @description Calculates the Concentration Index (CI) on a daily 
#' precipitation serie.
#' @usage ci(\dots, interval = 1)
#' @param \dots a set of daily precipitation series.
#' @param interval the interval in millimeters applied for calculating the 
#' concentration index. (Default value: 1)
#' @return A data.frame containing the following variables:
#' \itemize{
#' \item \code{dataset} is the precipitation serie name.
#' \item \code{a} is the constant of the exponential curve defined via the 
#' least-squares method.
#' \item \code{b} is the constant of the exponential curve defined via the 
#' least-squares method.
#' \item \code{r2} is the determination of coefficient of determination (R2) of 
#' the exponential model.
#' \item \code{A} is the area over the exponential curve.
#' \item \code{S} is the are compressed for the exponential curve.
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
#' ci(daily, interval = 1)
#' @keywords concentration index precipitation
#' @export
ci <- function(..., interval = 1) precintcon.ci.analysis(..., interval = interval, args = as.character(match.call()[1:length(list(...))+1]))