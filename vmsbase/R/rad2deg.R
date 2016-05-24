
#' Radians to Degree Angles conversion function
#' 
#' \code{rad2deg} converts radian heading data to degree values.
#' 
#' @param heading    The value in radians of the heading.
#'
#' @return The function returns the heading value converted in degrees
#' 
#' @usage rad2deg(heading)
#' 
#' @export rad2deg
#' 

#Heading edit
#Edits the heading data

rad2deg <- function (heading)
{
  
  heading <- heading * pi* (2/360)
  
  return(heading)
  
}