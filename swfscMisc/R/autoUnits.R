#' @title Auto Time Interval Units
#' @description Convert time interval units to natural values based on 
#'   magnitude of difference.
#' 
#' @param x an object inheriting from class \code{\link{difftime}}
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' autoUnits(as.difftime("0:3:35"))
#' autoUnits(as.difftime("15:3:35"))
#' autoUnits(ISOdate(2000, 5, 1) - ISOdate(2000, 4, 20))
#' 
#' @export
#' 
autoUnits <- function(x) {
  if(!inherits(x, "difftime")) stop("'x' is not a difftime object")
  if(is.na(x)) return(x)
  units(x) <- "secs"
  if(x > 60) {
    units(x) <- "mins"
    if(x > 60) {
      units(x) <- "hours"
      if(x > 24) {
        units(x) <- "days"
        if(x > 7) {
          units(x) <- "weeks"
        }
      }
    }
  }
  return(x)
}