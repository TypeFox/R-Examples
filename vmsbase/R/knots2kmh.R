#' Knots to Km/h speed conversion function
#' 
#' \code{kno2kmh} converts speed in knots to km/h values.
#' 
#' @param speed    The value in knots of the speed.
#'
#' @return The function returns the speed value converted in km/h
#' 
#' @usage kno2kmh(speed)
#'  
#' @export kno2kmh
#' 

#KNOTS2KMH
#Converts speed from knots to km/h

kno2kmh <- function (speed)
  {
  
    speed <- speed * 1.852
    
    return(speed)
  
  }