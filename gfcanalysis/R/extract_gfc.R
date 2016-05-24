#' @import methods
#' @import rgdal
#' @import raster
#' @importFrom sp bbox spTransform CRS proj4string
make_tile_mosaic <- function(aoi, data_folder, data_year, filename="",
                             stack="change", ...) {
    if (stack == 'change') {
        image_names <- c('treecover2000', 'loss', 'gain', 'lossyear', 
                         'datamask')
        band_names <- image_names
    } else if (stack == 'first') {
        image_names <- 'first'
        band_names <- c('Band3', 'Band4', 'Band5', 'Band7')
    } else if (stack == 'last') {
        image_names <- 'last'
        band_names <- c('Band3', 'Band4', 'Band5', 'Band7')
    } else {
        stop('"stack" must be equal to "change", "first", or "last"')
    }

    tiles <- calc_gfc_tiles(aoi)
    # Transform aoi to match tiles CRS so it can be used later for cropping
    aoi <- spTransform(aoi, CRS(proj4string(tiles)))
    file_root <- paste0('Hansen_GFC', data_year, '_')

    tile_stacks <- c()
    for (n in 1:length(tiles)) {
        tile <- tiles[n]
        min_x <- bbox(tile)[1, 1]
        max_y <- bbox(tile)[2, 2]
        if (min_x < 0) {
            min_x <- paste0(sprintf('%03i', abs(min_x)), 'W')
        } else {
            min_x <- paste0(sprintf('%03i', min_x), 'E')
        }
        if (max_y < 0) {
            max_y <- paste0(sprintf('%02i', abs(max_y)), 'S')
        } else {
            max_y <- paste0(sprintf('%02i', max_y), 'N')
        }
        file_suffix <- paste0('_', max_y, '_', min_x, '.tif')
        filenames <- file.path(data_folder, paste0(file_root, image_names, 
                                                   file_suffix))
        tile_stack <- crop(stack(filenames), aoi, datatype='INT1U', 
                           format='GTiff', options="COMPRESS=LZW")
        names(tile_stack) <- band_names
        tile_stacks <- c(tile_stacks, list(tile_stack))
    }

    if (length(tile_stacks) > 1) {
        # See http://bit.ly/1dJPIeF re issue in raster that necessitates below 
        # workaround TODO: Contact Hijmans re possible fix
        mosaic_list <- function(x, fun, datatype, format, options, overwrite, 
                                tolerance=0.05, filename="") {
            mosaic_args <- x
            if (!missing(fun)) mosaic_args$fun <- fun
            if (!missing(tolerance)) mosaic_args$tolerance <- tolerance
            if (!missing(datatype)) mosaic_args$datatype <- datatype
            if (!missing(format)) mosaic_args$format <- format
            if (!missing(options)) mosaic_args$options <- options
            if (!missing(overwrite)) mosaic_args$overwrite <- overwrite
            mosaic_args$filename <- filename
            do.call(mosaic, mosaic_args)
        }
        tile_mosaic <- mosaic_list(tile_stacks, fun='mean', filename=filename, 
                                   datatype='INT1U', format='GTiff', 
                                   options='COMPRESS=LZW', ...)
    } else {
        tile_mosaic <- tile_stacks[[1]]
        if (filename != '') {
            tile_mosaic <- writeRaster(tile_mosaic, filename=filename, 
                                       datatype="INT1U", format="GTiff", 
                                       options="COMPRESS=LZW", ...)
        }
    }
    names(tile_mosaic) <- band_names
    NAvalue(tile_mosaic) <- -1

    return(tile_mosaic)
}

#' Scale the first or last top of atmosphere (TOA) reflectance images
#'
#' This function applies the scale factors provided by Hansen et al. to rescale 
#' the first and last TOA reflectance images from integer to floating point.  
#' The following scale factors are used: band 3, 508; band 4, 254; band 5, 363; 
#' band 7, 423.  The output datatype is FLT4S.
#'
#' @export
#' @param x the "first" or "last" image for a given aoi as a \code{RasterStack} 
#' (see \code{stack} option for \code{\link{extract_gfc}}).
#' @param ... additional arguments as for \code{\link{writeRaster}}, such as 
#' \code{filename}, or \code{overwrite}.
#' @seealso \code{\link{extract_gfc}}
#' @return \code{RasterStack} of TOA reflectance values
scale_toar <- function(x, ...) {
    if (!nlayers(x) == 4) {
        stop('input image should have 4 bands')
    }
    scale_func <- function(b3, b4, b5, b7) {
        b3 <- (b3 - 1) / 508
        b4 <- (b4 - 1) / 254
        b5 <- (b5 - 1) / 363
        b7 <- (b7 - 1) / 423
        return(cbind(b3, b4, b5, b7))
    }
    x <- overlay(x, fun=scale_func, datatype='FLT4S',
                 format='GTiff', options="COMPRESS=LZW", ...)
    return(x)
}

#' Extracts GFC data for a given AOI
#'
#' This function extracts a dataset for a given AOI from a series of 
#' pre-downloaded GFC tiles. The \code{\link{download_tiles}} function should 
#' be used beforehand in order to download the necessary data to the specified
#' \code{data_folder}. Note that the output file format is fixed as GeoTIFF 
#' with LZW compression.
#'
#' The \code{stack} option can be "change" (the default), "first", or "last".  
#' When set to "change", the forest change layers (treecover2000, loss, gain, 
#' lossyear, and datamask) will be extracted for the given \code{aoi}. The 
#' "first" and "last" options will mosaic the 2000 or last year composite top 
#' of atmosphere (TOA) reflectance images (respectively).
#'
#' @seealso \code{\link{download_tiles}}, \code{\link{annual_stack}}, 
#' \code{\link{gfc_stats}}
#'
#' @export
#' @import rgdal
#' @importFrom sp spTransform CRS proj4string proj4string<-
#' @importFrom rgeos gBuffer
#' @param aoi an Area of Interest (AOI) as a \code{SpatialPolygons*} object.  
#' If the AOI is not in WGS 1984 (EPSG:4326), it will be reprojected to WGS84.
#' @param data_folder folder where downloaded GFC product tiles are located 
#' (see \code{\link{download_tiles}} function.
#' @param to_UTM if TRUE, then reproject the output into the UTM zone of the 
#' AOI centroid. If FALSE, retain the original WGS84 projection of the GFC 
#' tiles.
#' @param stack the layers to extract from the GFC product. Defaults to 
#' "change". See Details.
#' @param data_year which version of the Hansen data to use
#' @param ... additional arguments as for \code{\link{writeRaster}}, such as 
#' \code{filename}, or \code{overwrite}.
#' @return \code{RasterStack} with GFC layers
extract_gfc <- function(aoi, data_folder, to_UTM=FALSE, stack="change", 
                        data_year=2015, ...) {
    if (stack == 'change') {
        band_names <- c('treecover2000', 'loss', 'gain', 'lossyear', 
                        'datamask')
    } else if (stack == 'first') {
        band_names <- c('Band3', 'Band4', 'Band5', 'Band7')
    } else if (stack == 'last') {
        band_names <- c('Band3', 'Band4', 'Band5', 'Band7')
    } else {
        stop('"stack" must be equal to "change", "first", or "last"')
    }

    if (to_UTM) {
        tile_mosaic <- make_tile_mosaic(aoi, data_folder, stack=stack, 
                                        data_year=data_year, ...)
        # Project to UTM for plotting and analysis of change in forest area.  
        # Calculate UTM zone based on bounding polygon of tile mosaic.
        bounding_poly <- as(extent(tile_mosaic), "SpatialPolygons")
        proj4string(bounding_poly) <- proj4string(tile_mosaic)
        utm_proj4string <- utm_zone(bounding_poly, proj4string=TRUE)
        # Use nearest neighbor since the data is categorical
        band_name <- names(tile_mosaic)
        tile_mosaic <- projectRaster(tile_mosaic, crs=utm_proj4string, 
                                     method='ngb', datatype='INT1U', 
                                     format='GTiff', options="COMPRESS=LZW", 
                                     ...)
        names(tile_mosaic) <- band_names
        NAvalue(tile_mosaic) <- -1
    } else {
        tile_mosaic <- make_tile_mosaic(aoi, data_folder, stack=stack, 
                                        data_year=data_year, ...)
        NAvalue(tile_mosaic) <- -1
    }

    return(tile_mosaic)
}
