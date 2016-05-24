#' @include precintcon.stat.analysis.r
NULL

#' @name stat
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com} 
#' @aliases precintcon.stat.analysis stat 
#' @title Basic statistics for precipitation datasets. 
#' @description Calculates of basic statistics of precipitation datasets. 
#' @usage stat(\dots) 
#' @param \dots a set of daily or monthly precipitation serie.
#' @return A data.frame (precintcon.stat) containing the following variables:
#' \itemize{ 
#' \item \code{dataset} is the precipitation serie name.
#' \item \code{mean.daily} is the daily average precipitation. 
#' It exists only for daily datasets.
#' \item \code{sd.daily} is the standard deviation of a daily precipitation serie. 
#' It exists only for daily datasets.
#' \item \code{var.daily} is the variance of a daily precipitation serie. 
#' It exists only for daily datasets.
#' \item \code{mean.monthly} is the monthly average precipitation.
#' \item \code{sd.monthly} is the standard deviation of a monthly precipitation serie.
#' \item \code{var.monthly} is the variance of a monthly precipitation serie.
#' \item \code{total} is the total precipitation.
#' }
#' @examples 
#' ##
#' # Loading the monthly precipitation serie.
#' data(monthly)
#' 
#' ##
#' # Plotting the deciles.
#' stat(monthly)
#' @keywords summary precipitation 
#' @export
stat <- function(...) precintcon.stat.analysis(..., args = as.character(match.call()[1:length(list(...))+1]))