#' Create grid of locations.
#' 
#' \code{create.pgrid2} creates a grid of locations fusing vectors of x and y coordinates.
#' 
#' The key argument in the function midpoints. If this is \code{TRUE}, it is assumed that the boundaries of the spatial domain correspond to the midpoints of the cell/pixel in the grid. Otherwise, it is assumed that the boundaries correspond to the actual borders of the region of interest. If \code{poly.coords} is supplied, the grid returned is the grid of midpoints contained in the convex hull of \code{poly.coords}.
#' 
#' @param xgrid A vector of locations in the x direction.
#' @param ygrid A vector of location in the y direction.
#' @param midpoints A logical value (\code{TRUE} or \code{FALSE}) indicating whether the boundary values are for the midpoint of a pixel (\code{midpoints = TRUE}) or for the boundary of the spatial domain in general (\code{midpoints = FALSE}, in which case the midpoints are calculated internally). Default is \code{FALSE}.
#' @param poly.coords An \eqn{n \times 2} matrix with the coordinates specifying the polygon vertices of the true spatial domain of interest within the rectangular boundaries provided by \code{xmin}, \code{xmax}, \code{ymin}, and \code{ymax}. If this is provided, the \code{pgrid} returned will be within the convex hull of \code{poly.coords}.
#' 
#' @return Returns an object of class pgrid with the following components: 
#' \item{pgrid}{An \eqn{n \times 2} matrix of locations (the midpoints of the pixelized grid).}
#' \item{m}{The number of rows in pgrid.}
#' \item{p.in.grid}{A vector of 0s and 1s indicating whether the midpoint of each pixel is in the convex hull of \code{poly.coords}. If \code{poly.coords} is not provided, this is a vector of 1s.}
#' \item{ubx}{The pixel boundaries in the x-direction.}
#' \item{uby}{The pixel boundaries in the y-direction.}
#' \item{upx}{The pixel midpoints in the x-direction.}
#' \item{upy}{The pixel midpoints in the y-direction.}
#' 
#' @author Joshua French
#' @importFrom splancs inout
#' @export
#' @examples 
#' seq1 = seq(0, 1, len = 101)
#' pgrida <- create.pgrid2(seq1, seq1, midpoint = FALSE)
#' seq2 = seq(.005, .995, len = 100)
#' pgridb <- create.pgrid2(seq2, seq2, midpoint = TRUE)
#' # pgrids produced match
#' range(pgrida$pgrid - pgridb$pgrid)
create.pgrid2 <- function(xgrid, ygrid, midpoints = FALSE, poly.coords = NULL)
{
  nx <- length(xgrid)
  ny <- length(ygrid)
  
  if(!midpoints)
  {
    #Rename xgrid and ygrid (for later consistency).  
    ubx <- xgrid; uby <- ygrid
    
    #Determine x and y midpoints of all of the pixels
    upx <- (xgrid[-nx] + xgrid[-1])/2
    upy <- (ygrid[-ny] + ygrid[-1])/2
  }else
  {
    #Rename xgrid and ygrid (for later consistency)
    upx <- xgrid; upy <- ygrid
    
    #Create boundaries for pixels
    ubx <- (xgrid[-nx] + xgrid[-1])/2
    ubx <- c(2 * upx[1] - ubx[1], ubx, 2 * upx[nx] - ubx[nx - 1])
    
    uby <- (ygrid[-ny] + ygrid[-1])/2
    uby <- c(2 * upy[1] - uby[1], uby, 2 * upy[ny] - uby[ny - 1])
  }
  
  #If coords are supplied, create pgrid based on whether points
  #are contained in polygon of poly.coords.
  if(!is.null(poly.coords))
  {
    all.grid <- as.matrix(expand.grid(upx, upy))
    #Determine points of rectangular grid (based on xgrid and ygrid)
    #within poly.coords
    pip <- inout(all.grid, poly.coords, bound = TRUE)
    #Extract prediction coordinates within border 
    pgrid <- as.matrix(all.grid[pip == 1,])
    #Determine number of prediction locations
    np <- nrow(pgrid)
    #Determine which points are a prediction coordinates within 
    #rectangular grid
    p.in.grid <- (pip == 1)
  }else
  {
    pgrid <- as.matrix(expand.grid(upx, upy))
    np <- length(upx) * length(upy)
    p.in.grid <- rep(TRUE, np)
  }
  
  out <- (list(pgrid = as.matrix(pgrid), np = np, p.in.grid = p.in.grid, 
               ubx = ubx, uby = uby, upx = upx, upy = upy))	
  class(out) <- "pgrid"
  return(out)
}
