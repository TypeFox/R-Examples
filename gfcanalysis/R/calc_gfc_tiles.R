#' Calculate the GFC product tiles needed for a given AOI
#'
#' Intersects an (optionally buffered) AOI with the GFC 
#' product grid to determine what tiles are need to cover the AOI.
#'
#' @export
#' @import rgdal
#' @importFrom sp spTransform CRS proj4string
#' @importFrom rgeos gBuffer gIntersects gConvexHull
#' @param aoi an Area of Interest (AOI) as a \code{SpatialPolygons*} object.  
#' If the AOI is not in the WGS84 geographic coordinate system, it will be 
#' reprojected to WGS84.
#' @return a \code{SpatialPolygonsDataFrame} of the GFC tiles needed to cover 
#' the AOI
#' @examples
#' tiles <- calc_gfc_tiles(test_poly)
#' plot(tiles)
#' plot(test_poly, lt=2, add=TRUE)
calc_gfc_tiles <- function(aoi) {
    aoi <- spTransform(aoi, CRS(proj4string(gfc_tiles)))
    intersecting <- as.logical(gIntersects(gfc_tiles, 
                                           gConvexHull(aoi), 
                                           byid=TRUE))
    if (sum(intersecting) == 0) {
        stop('no intersecting GFC tiles found')
    } else {
        gfc_tiles <- gfc_tiles[intersecting, ]
    }
    return(gfc_tiles)
}
