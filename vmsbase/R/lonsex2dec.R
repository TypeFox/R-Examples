#' Sexagesimal to Decimal Longitude conversion function
#' 
#' \code{lonsex2dec} converts sexagesimal longitude data to decimal longitude values.
#' 
#' @param degree    The value in degrees of the longitude.
#' @param minute    The value in minutes of the longitude.
#' @param second    The value in seconds of the longitude.
#' @param direction   The direction, as "E" or "W", of the longitude.
#'
#' @return The function returns the longitude value converted in the decimal numeral system
#' 
#' @usage lonsex2dec(degree, minute, second, direction)
#' 
#' @export lonsex2dec
#' 

#LON_SEX2DEC
#Converts longitude from sexagesimal to decimal numeral system.

lonsex2dec <- function (degree, minute, second, direction)
{
  if(direction == "E")
  {
    declon <- degree+(minute/60)+(second/3600)
  }else{
    declon <- -(degree+(minute/60)+(second/3600))
  }
  return(declon)
}