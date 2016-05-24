#' @name ff.index
#' @aliases ff.index
#' @author Lucas Venezian Povoa
#' 
#' @title F factor
#' @description Calculates the Approximated Rainfall Erosivity Factor according to the ff index.
#' @usage ff.index(object)
#' @param object is a daily or monthly precipitation serie
#' @return the ff index in millimeters
#' @examples 
#' ##
#' # Loading the daily precipitation serie.
#' data(daily)
#' 
#' ##
#' # Calculating ff index
#' ff.index(daily)
#' @references
#' Ferro, V., Giordano, G., Iovino, M. (1991). Isoerosivity and Erosion Risk Map for Sicily.
#' Hydrolog. Sci. J. 36, 549-564
#' @export
ff.index <- function(object) {
  
  m <- as.monthly(object)
  
  y <- as.annual(object)
  
  result <- 0;
  
  for (j in 1:nrow(y)) {
    my <- m[m$year == y$year[j],]
    for (i in 1:nrow(my))
      result <- result + ((my$precipitation[i]**2) / y$precipitation[j])
  }
  
  return(result / nrow(y))
}