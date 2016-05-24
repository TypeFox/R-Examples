#' @title Get configuration of a light cone (LC)
#'
#' @description 
#' \code{get_LC_config} obtains the PLC or FLC at a particular 
#' \eqn{(\mathbf{r}, t)} from a \eqn{(N+1)D} field based on the LC 
#' template from \code{\link{compute_LC_coordinates}} (or 
#' \code{\link{setup_LC_geometry}}).
#' 
#' @param coord space-time coordinate \eqn{(\mathbf{r}, t)}
#' @param field spatio-temporal field; either a matrix or a 3-dimensional array 
#' with time \eqn{t} as the first coord, and the spatial coords in order.  
#' Make sure to see also \code{\link{compute_LC_coordinates}} for correct 
#' formatting.
#' @param LC.coordinates template coords for the LC
#' @keywords method utilities
#' @export
#' @seealso \code{\link{compute_LC_coordinates}}
#' @examples
#' AA = matrix(rnorm(40), ncol = 5)
#' image2(AA)
#' LCind = compute_LC_coordinates(speed = 1, horizon = 1, shape = "cone")
#' AA
#' get_LC_config(cbind(5,2), AA, LCind)
#' # a time series example
#' data(nhtemp)
#' xx <- c(nhtemp)
#' LCind = compute_LC_coordinates(speed = 1, horizon = 4, shape = "cone", space.dim = 0)
#' cc <- get_LC_config(6, xx, LCind)
#' 

get_LC_config <- function(coord, field, LC.coordinates) {
  
  coord <- cbind(coord)
  if (length(dim(coord)) == 0) {
    #stopifnot(length(coord) == ncol(LC.coordinates))
    if (ncol(LC.coordinates) == 1) {
      result <- field[coord + LC.coordinates]
    } else { # if greater than one, use sweep
      result <- field[sweep(LC.coordinates, 2, coord, "+")]
    }
  } else {
    if (ncol(LC.coordinates) == 1) {
      result <- 
        sapply(coord, 
               function(x) {
                 return(field[x + LC.coordinates])
               })
    } else {
      result <- apply(coord, 1, 
                      function(x) {
                        return(field[sweep(LC.coordinates, 2, x, "+")])
                      })
    }
  }
  return(result)
}

