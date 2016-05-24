#' Calculates the pixel area for each line of a raster
#'
#' @import raster
#' @importFrom geosphere areaPolygon
#' @param x a \code{Raster*} object
#' @return a vector with the area in square meters of the pixels in each line 
#' of \code{x} (vector has length equal to \code{nrow(x)})
calc_pixel_areas <- function(x) {
    # Construct polygons for a single column of the raster
    xleft <- xmin(x)
    xright <- xmin(x) + xres(x)
    # Note that ylower and yupper are setup so that the polygons are ordered
    # from high latitude to lowest latitude - needed because rasters are 
    # addressed starting from the upper left corner (highest latitude)
    ylower <- seq(from=(ymax(x) - yres(x)), by=(-yres(x)), length.out=nrow(x))
    yupper <- seq(from=ymax(x), by=(-yres(x)), length.out=nrow(x))
    poly_areas <- function(xl, xr, yl, yu) {
        areaPolygon(matrix(c(xl, yl,
                             xr, yl,
                             xr, yu,
                             xl, yu), ncol=2, byrow=TRUE))
    }
    pixel_areas <- mapply(poly_areas, xleft, xright, ylower, yupper)
}
    
#' Scales a raster by the area of each pixel in meters
#'
#' Calculates the area (in meters) of each pixel in a raster, scales the value 
#' of each pixel by the area, applies the desired scale factor, and returns the 
#' result as a \code{RasterLayer}.  Useful for calculating class areas based on 
#' a classified raster in a geographic coordinate system.  Assumes that raster 
#' is not rotated (latitudes of every pixel in a given row are identical).  
#' Processes block by block to support handling very large rasters.
#'
#' @import raster
#' @param x a \code{Raster*} object
#' @param pixel_areas a vector giving the area of each cell in a single column 
#' of the raster. See \code{\link{calc_pixel_areas}}. If NULL, this vector will 
#' be calculated based on the coordinate system of \code{x}.
#' @param scale_factor a value to scale the results by
#' @param filename (optional) filename for output raster
#' @param datatype (optional) datatype for output raster see 
#' \code{\link{dataType}} NOT YET SUPPORTED
#' @return \code{RasterLayer} with pixel areas (in meters)
scale_by_pixel_area <- function(x, filename, datatype, pixel_areas=NULL,
                                scale_factor=1) {
    if (is.null(pixel_areas)) {
        pixel_areas <- calc_pixel_areas(x)
    }
    if (class(x) == 'RasterLayer') {
        out <- raster(x)
    } else if (class(x) %in% c('RasterStack', 'RasterBrick')) {
        out <- brick(x, values=FALSE)
    } else {
        stop('x must be a Raster* object')
    }
    # Write areas by block to avoid running out of memory with very large 
    # rasters
    bs <- blockSize(out)
    if (missing(filename)) filename <- rasterTmpFile()
    out <- writeStart(out, filename=filename)
    for (block_num in 1:bs$n) {
        start_row <- bs$row[block_num]
        nrows <- bs$nrows[block_num]
        x_bl <- getValuesBlock(x, start_row, nrows)
        these_pixel_areas <- rep(pixel_areas[start_row:(start_row + nrows - 1)], each=ncol(x))
        if (is.matrix(x_bl)) {
            # Block of a stack or brick has more than one column
            x_bl <- apply(x_bl, 2, function(col_vals) col_vals * 
                          these_pixel_areas * scale_factor)
        } else {
            # Block of a rasters is just a vector
            x_bl <- x_bl * these_pixel_areas * scale_factor
        }
        writeValues(out, x_bl, start_row)
    }
    out <- writeStop(out)
    return(out)
}
