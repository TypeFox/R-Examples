#' @include precintcon.spi.per.year.analysis.r
NULL

#' @name spi.per.year 
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com}
#' @aliases precintcon.spi.per.year.analysis spi.per.year 
#' @title Standardized Precipitation Index 
#' @description Calculates the Standardized Precipitation Index (SPI) per 
#' year on a daily or monthly precipitation serie. 
#' @usage spi.per.year(object, period = 3, distribution = "Gamma", FUN = mean) 
#' @param object a daily or monthly precipitation serie.
#' @param period the number of months to be aggregate in the calculation 
#'                 of the standardized precipitation index. (Default value: 3)
#' @param distribution (it has no effect yet).
#' @param FUN the function used to summarize the standardized 
#'              precipitation index per year. (Default function: mean).
#' @return A data.frame (precintcon.spi.per.year) containing the following variables:
#' \itemize{
#' \item \code{year} is the year.
#' \item \code{spi} is the standardized precipitation index.
#' }
#' @seealso
#' \code{\link{spi}}
#' \code{\link{read.data}}
#' @examples 
#' ##
#' # Loading the monthly precipitation serie.
#' data(monthly)
#' 
#' ##
#' # Performing the Standardized Precipitation Index
#' spi.per.year(monthly, period = 3)
#' @keywords rainfall anomaly index precipitation 
#' @export
spi.per.year <- function(object, period = 3, distribution = "Gamma", FUN = mean) 
  precintcon.spi.per.year.analysis(object, period, distribution, FUN)