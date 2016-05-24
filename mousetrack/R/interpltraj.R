#' Interpolate trajectories to a fix number of time bins
#' 
#' @param x x coordinate
#' @param y y coordinate.
#' @param singlepoint a logical flag to indicate whether interpolation is done on a single coordinate point (TRUE) or two points (FALSE)
#' @param tsmax the new length of the interpolated trajectory.
#' @return The interpolated trajectory.
#' @examples
#'
#' data(mousetrack)
#' x = set$x; y = set$y;
#' singlepoint = FALSE; tsmax = 101
#' interpltraj(x, y, singlepoint, tsmax)

.packageName <- 'mousetrack'

interpltraj <- function (x, y, singlepoint, tsmax){
    
    if (singlepoint == TRUE) {
            
        interts = approx(x, y = NULL, n = tsmax)$y            
            
    } else { ## cases with x-y coordinates

        ans = approx(x, y, n = tsmax)                
        interts =  cbind(ans$x, ans$y, deparse.level = 0)
                
    }
    
    return(interts)

}
