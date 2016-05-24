#' Create grid of locations.
#' 
#' \code{create.pgrid} creates a grid of locations from the boundaries of domain and other information.
#' 
#' The key argument in the function midpoints. If this is \code{TRUE}, it is assumed that the boundaries of the spatial domain correspond to the midpoints of the cell/pixel in the grid. Otherwise, it is assumed that the boundaries correspond to the actual borders of the region of interest. If \code{poly.coords} is supplied, the grid returned is the grid of midpoints contained in the convex hull of \code{poly.coords}.
#' 
#' @param xmin The minimum value of the boundary of the x coordinates of the spatial domain.
#' @param xmax The maximum value of the boundary of the x coordinates of the spatial domain.
#' @param ymin The minimum value of the boundary of the y coordinates of the spatial domain.
#' @param ymax The maximum value of the boundary of the y coordinates of the spatial domain.
#' @param nx The number of gridpoints/cells/pixels in the x direction.
#' @param ny The number of gridpoints/cells/pixels in the y direction.
#' @param midpoints A logical value (\code{TRUE} or \code{FALSE}) indicating whether the boundary values are for the midpoint of a pixel (\code{midpoints = TRUE}) or for the boundary of the spatial domain in general (\code{midpoints = FALSE}), in which case the midpoints are calculated internally). Default is \code{FALSE}.
#' @param poly.coords An \eqn{n \times 2} matrix with the coordinates specifying the polygon vertices of the true spatial domain of interest within the rectangular boundaries provided by \code{xmin}, \code{xmax}, \code{ymin}, and \code{ymax}. If this is provided, the \code{pgrid} returned will be within the convex hull of \code{poly.coords}.
#' 
#' @return Returns an object of class \code{pgrid} with the following components: 
#' \item{pgrid}{An \eqn{n \times 2} matrix of locations (the midpoints of the pixelized grid).}
#' \item{m}{The number of rows in pgrid.}
#' \item{p.in.grid}{A vector of 0s and 1s indicating whether the midpoint of each pixel is in the convex hull of \code{poly.coords}. If \code{poly.coords} is not provided, this is a vector of 1s.}
#' \item{ubx}{The pixel boundaries in the x direction.}
#' \item{uby}{The pixel boundaries in the y direction.}
#' \item{upx}{The pixel midpoints in the x direction.}
#' \item{upy}{The pixel midpoints in the y direction.}
#' 
#' @author Joshua French
#' @importFrom splancs inout
#' @export
#' @examples 
#' pgrida <- create.pgrid(0, 1, 0, 1, nx = 50, ny = 50, midpoints = FALSE)
#' pgridb <- create.pgrid(.01, .99, .01, .99, nx = 50, ny = 50, midpoints = TRUE)
create.pgrid <- function(xmin, xmax, ymin, ymax, nx, ny, midpoints = FALSE,
                         poly.coords = NULL)
{
  if(midpoints)
  {
    
    xstep <- (xmax-xmin)/(nx - 1)	#Calculate the pixel width
    ystep <- (ymax-ymin)/(ny - 1)	#Calculate the pixel height
    
    #Determine x and y midpoints of all of the pixels	
    upx <- xmin + 0:(nx - 1) * xstep 
    upy <- ymin + 0:(ny - 1) * ystep 
    
    #Create boundaries for pixels
    ubx <- xmin + 0:nx * xstep - xstep/2
    uby <- ymin + 0:ny * ystep - ystep/2
  }
  else
  {
    xstep <- (xmax-xmin)/nx	#Calculate the pixel width
    ystep <- (ymax-ymin)/ny	#Calculate the pixel height
    
    #Determine x and y midpoints of all of the pixels	
    upx <- xmin + xstep/2 + 0:(nx-1) * xstep 
    upy <- ymin + ystep/2 + 0:(ny-1) * ystep 
    
    #Create boundaries for pixels
    ubx <- xmin + 0:nx * xstep
    uby <- ymin + 0:ny * ystep
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
