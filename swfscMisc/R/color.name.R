#' @title Color Name
#' @description Return the name of a color listed given the number.
#'  
#' @param i integer specifying color .
#' 
#' @return character value of 'i' color.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @importFrom grDevices colors
#' @export
#' 
color.name <- function(i) colors()[i]