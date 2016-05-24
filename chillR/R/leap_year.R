#' Leap year finder
#' 
#' This function determines whether a given year is a leap year
#' 
#' Takes a year number as input, and returns TRUE if this is a leap year, and
#' FALSE if not
#' 
#' @param x integer value, representing year number
#' @return boolean variable (TRUE or FALSE)
#' @author Eike Luedeling, but based on pseudocode from Wikipedia
#' @references https://en.wikipedia.org/wiki/Leap_year
#' @keywords utility
#' @examples
#' 
#' 
#' leap_year(2015)  
#' leap_year(2016)
#' 
#'  
#' @export leap_year
leap_year<-function(x) if(!x/4==trunc(x/4)) FALSE else
                         if(!x/100==trunc(x/100)) TRUE else
                           if(!x/400==trunc(x/400)) FALSE else TRUE
