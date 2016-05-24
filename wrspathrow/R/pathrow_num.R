load_wrs_data <- function(wrs_type, wrs_mode) {
    if (wrs_type == 2) {
        wrs_polys <- wrs2_asc_desc
    } else if (wrs_type == 1) {
        wrs_polys <- wrs1_asc_desc
    } else {
        stop('wrs_type must be 1 or 2')
    }
    if (!(wrs_mode %in% c('D', 'A'))) {
        stop('wrs_mode must be "D", "A" or c("D", "A")')
    }
    return(wrs_polys[wrs_polys@data$MODE %in% wrs_mode, ])
}

intersect_wrs_polys <- function(wrs_polys, x, as_polys) {
    intersecting <- as.logical(gIntersects(wrs_polys, x, byid=TRUE))
    if (sum(intersecting) == 0) {
        stop('no intersecting pathrows found')
    } else {
        wrs_polys <- wrs_polys[intersecting, ]
        wrs_polys <- wrs_polys[order(wrs_polys$PATH, wrs_polys$ROW), ]
        if (!as_polys) {
            wrs_polys <- data.frame(PATH=wrs_polys@data$PATH, ROW=wrs_polys@data$ROW)
        }
        return(wrs_polys)
    }
}

#' Get WRS-2 path/row numbers for a given spatial object
#'
#' @export
#' @docType methods
#' @import methods
#' @import wrspathrowData
#' @rdname pathrow_num-methods
#' @param x a spatial object
#' @param wrs_type 1 (for WRS-1) or 2 (for WRS-2)
#' @param wrs_mode either 'D' for descending (daytime) or 'A' for ascending 
#' @param as_polys if FALSE (default) return a data.frame. If TRUE, return a 
#' \code{SpatialPolygonsDataFrame}.
#' @return data.frame with path and row as integers, or, if as_polys=TRUE, a 
#' \code{SpatialPolygonsDataFrame}
#' @examples
#' \dontrun{
#' library(sp)
#'
#' pathrow_num(test_poly)
#'
#' x <- pathrow_num(test_poly, as_polys=TRUE)
#' plot(x)
#' plot(test_poly, add=TRUE, lty=2, col="#00ff0050")
#' text(coordinates(x), labels=paste(x$PATH, x$ROW, sep=', '))
#' }
setGeneric("pathrow_num", function(x, wrs_type='2', wrs_mode='D', 
                                   as_polys=FALSE) {
    standardGeneric("pathrow_num")
})

#' @rdname pathrow_num-methods
#' @importFrom raster extent projectExtent crs
#' @importFrom rgeos gIntersects
#' @aliases pathrow_num,Raster-method
setMethod("pathrow_num", signature(x="Raster"),
    function(x, wrs_type, wrs_mode, as_polys) {
        wrs_polys <- load_wrs_data(wrs_type, wrs_mode)
        x_wgs84 <- projectExtent(x, crs=crs(wrs_polys))
        x_wgs84_sp <- as(extent(x_wgs84), 'SpatialPolygons')
        return(intersect_wrs_polys(wrs_polys, x_wgs84_sp, as_polys))
    }
)

#' @rdname pathrow_num-methods
#' @importFrom rgeos gIntersects
#' @importFrom sp CRS proj4string spTransform
#' @import rgdal
#' @aliases pathrow_num,Spatial-method
setMethod("pathrow_num", signature(x="Spatial"),
    function(x, wrs_type, wrs_mode, as_polys) {
        wrs_polys <- load_wrs_data(wrs_type, wrs_mode)
        x_wgs84 <- spTransform(x, CRS(proj4string(wrs_polys)))
        return(intersect_wrs_polys(wrs_polys, x_wgs84, as_polys))
    }
)
