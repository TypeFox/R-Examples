#' Generate a spatial grid
#'
#' Produces an arbitrary grid in any user-defined coordinate system. Used by 
#' \code{gfcanalysis} for producing the 10x10 degree WGS84 grid that the GFC 
#' product is tiled on, so that an AOI polygon can be intersected with the grid 
#' to calculate the appropriate tiles to download.
#'
#' @seealso \code{\link{download_tiles}}
#'
#' @export
#' @importFrom sp CRS GridTopology SpatialGrid
#' @param origin_x x coordinate of the origin
#' @param dx cell size in the x direction
#' @param max_x maximum value in x direction
#' @param origin_y y coordinate of the origin
#' @param dy cell size in the y direction
#' @param max_y maximum value in y direction
#' @param grid_proj4string coordinate system as a crs object (defaults to 
#' WGS-84)
#' @examples
#' gfc_tiles <- gen_grid(-180, 10, 180, -60, 10, 80)
gen_grid <- function(origin_x, dx, max_x, origin_y, dy, max_y,
                     grid_proj4string=NULL) {
    if (is.null(grid_proj4string)) grid_proj4string <- '+init=epsg:4326'
    # Based on code at http://bit.ly/1lfUOnV
    cells.dim <- c((max_x - origin_x) / dx,
                   (max_y - origin_y) / dy)
    gt <- GridTopology(c(origin_x+dx/2, origin_y+dy/2), c(dx, dy), cells.dim)
    grd <- SpatialGrid(gt, proj4string=grid_proj4string)
    spix <- as(grd, "SpatialPixels")
    spol <- as(spix, "SpatialPolygons")
    return(spol)
}
