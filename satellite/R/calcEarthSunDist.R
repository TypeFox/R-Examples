#' Compute earth-sun distance based on day of the year
#'
#' @description
#' The earth-sun distance for a particular day of the year is computed based on
#' one of several empirical formulas.
#' 
#' @param date Date of the sensor overpass; either a character string in a 
#' native date format (e.g. "YYYY-MM-DD", see \code{\link{as.Date}}) or a POSIX* 
#' object (see \code{\link{as.POSIXct}}). 
#' @param formula Formula to be applied, specified through the name of the 
#' author, i.e. one of "Spencer", "Mather" or "ESA".
#'
#' @return Numeric earth-sun distance (in AU).
#'
#' @export calcEarthSunDist
#' 
#' @details Computation of earth-sun distance using formulas provided by
#' Spencer (1971), Mather (2005) or ESA.
#' 
#' @references The formulas are taken from the following sources:
#' \itemize{
#'   \item Spencer: Spencer JW (1971) Fourier series representation of the position of 
#'   the sun. Search 2/5. Taken from 
#'   \url{https://goo.gl/lhi9UI}. 
#' 
#'   \item Mather: Paul M. Mather (2005) Computer Processing of Remotely-Sensed Images:
#'   An Introduction. Wiley, ISBN: 978-0-470-02101-9, 
#'   \url{http://eu.wiley.com/WileyCDA/WileyTitle/productCd-0470021012.html}.
#' 
#'   \item ESA: ESA Earth Observation Quality Control: Landsat frequently asked questions. 
#' }
#' See also: Bird R, Riordan C (1984) Simple solar spectral model for direct and 
#' diffuse irradiance on horizontal and tilted planes at the Earth's surface for
#' cloudless atmospheres. \url{http://www.nrel.gov/docs/legosti/old/2436.pdf}.
#' 
#' @examples
#' calcEarthSunDist(date = "2015-01-01", formula = "Spencer")
#' 
calcEarthSunDist <- function(date, formula = c("Spencer", "Mather", "ESA")){
  
  ## if not supplied, formula defaults to "Spencer"
  formula <- formula[1]
  
  day <- as.numeric(strftime(date, format = "%j"))
  if(formula == "Spencer"){
    pos <- 2 * pi * (day - 1) / 365
    (1/(1.000110 + 0.034221 * cos(pos) + 0.001280 * sin(pos) + 0.000719 * 
      cos(2 * pos) + 0.000077 * sin(2 * pos)))**0.5
  } else if(formula == "Mather"){
    1/(1 - 0.016729 * cos(0.9856 * (day - 4)))
  } else if(formula == "ESA"){
    1 - 0.016729 * cos((2 * pi) * (0.9856 * (day - 4) / 360))
  }
}

