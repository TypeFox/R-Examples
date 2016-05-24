#' Define a geographic origin
#'
#' \command{kissmigOrigin} is a support function to define the geographic origin for a \command{kissmig} call.
#' @usage
#' kissmigOrigin(grd, x, y, size)
#' @param grd a single RasterLayer as reference layer
#' @param x lower left x-coordinate of geographic origin
#' @param y lower left y-coordinate of geographic origin
#' @param size size of the quadratic origin
#'
#' @details
#' \command{kissmigOrigin} returns a rasterLayer characterized by the reference 
#' \command{grd}. Cell values are set to zero, except for cells of the origin defined by  
#' \command{extent(x, x+size, y, y+size)} which are set to one. 
#' @seealso \code{\link{kissmig}}
#' @export kissmigOrigin  

kissmigOrigin <- function(grd, x, y, size) {
  ans <- grd
  values(ans) <- 0
  values(ans)[cellsFromExtent(ans, extent(x, x+size, y, y+size))] <- 1.0
  return(ans)
}
