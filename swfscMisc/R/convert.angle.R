#' @title Angle Conversion
#' @description Converts angles between radians and degrees.
#' 
#' @param x numeric. The angle to be converted.
#' @param from,to character. Units to convert from and to. Can be 
#'   "radians" or "degrees" or any partial match (case-sensitive).
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' convert.angle(45, "deg", "rad")
#' convert.angle(4.5, "r", "d")
#' 
#' @export
#' 
convert.angle <- function(x, from = c("degrees", "radians"), to = c("radians", "degrees")) {
  from <- match.arg(from)
  to <- match.arg(to)
  
  x <- switch(from,
    degrees = x,
    radians = x * 180 / pi
  )
  
  result <- switch(to,
    degrees = x,
    radians = x * pi / 180
  )
  
  as.numeric(result)
}