#' Compute a geographic distance matrix
#' 
#' @author Bruno Vilela & Fabricio Villalobos
#' 
#' @description Calculates a geographic distance matrix based on a \code{PresenceAbsence} or 
#' a two column \code{matrix} of x(longitude) and y(latitude).
#' 
#' @param xy A \code{\link{PresenceAbsence}} object or a \code{matrix} with two columns (longitude, latitude).
#' @param asdist Logical, if \code{TRUE} the result will be an object of class \code{dist},
#' if \code{FALSE} the result will be an object of class \code{matrix}.
#' @param ... Arguments to be passed to the function \code{rdist.earth} of package fields.
#'   
#' @details This function basically facilitates the use of \code{rdist.earth} 
#' on a \code{PresenceAbsence} object, allowing also the user to have directly 
#' a \code{dist} object.
#' 
#' @return The user can choose between \code{dist} and \code{matrix} class object to be returned.
#' The resulting values are in kilometers (but see the argument 'miles' in \code{rdist.earth}).  
#'    
#' @examples \dontrun{
#' data(PAM)
#' distPAM <- lets.distmat(PAM)   
#' }
#' 
#' @export


lets.distmat <- function(xy, asdist = TRUE, ...) {
  
  # If a presence absence change to coordinates
  if (class(xy) == "PresenceAbsence"){
    xy <- xy[[1]][, 1:2]
  }
  
  # Default in Km
  if (exists("miles")) {
    # Calculate the distance matrix
    distan <- rdist.earth(xy, miles = FALSE, ...)
  } else {
    distan <- rdist.earth(xy, ...)
  }
  
  # Assuring that the diagonal of the matrix is zero
  diag(distan) <- 0
  
  # Transform in a distance matrix
  if (asdist) {
    distan <- as.dist(distan)
  }
  
  return(distan)
}

