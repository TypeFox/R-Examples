#' @include precintcon.spi.analysis.r
NULL

#' @name spi
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com} 
#' @aliases precintcon.spi.analysis spi 
#' @title Standardized Precipitation Index 
#' @description Calculates the Standardized Precipitation Index (SPI)
#' on a daily or monthly precipitation serie. 
#' @usage spi(object, period = 3, distribution = "Gamma") 
#' @param object a daily or monthly precipitation serie.
#' @param period the number of months to be aggregated in the calculation 
#' of the standardized precipitation index. (Default value: 3)
#' @param distribution it has no effect yet. (Default value: "Gamma")
#' @return A data.frame (precintcon.spi) containing the following variables:
#' \itemize{ 
#' \item \code{year} is the year.
#' \item \code{month} is the month.
#' \item \code{spi} is the standardized precipitation index.
#' }
#' @seealso 
#' \code{\link{precintcon.plot.spi}}
#' \code{\link{read.data}}
#' @examples 
#' ##
#' # Loading the monthly precipitation serie.
#' data(monthly)
#' 
#' ##
#' # Performing the Standardized Precipitation Index
#' spi(monthly, period = 3)
#' @keywords rainfall anomaly index precipitation
#' @export 
spi <- function(object, period = 3, distribution = "Gamma") precintcon.spi.analysis(object, period, distribution)