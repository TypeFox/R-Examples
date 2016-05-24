gen_year_list <- function(data_year) {
    if (data_year == 2013) {
        years <- seq(2000, 2012, 1)
    } else if (data_year == 2014) {
        years <- seq(2000, 2013, 1)
    } else if (data_year == 2015) {
        years <- seq(2000, 2014, 1)
    } else if (data_year > 2015) {
        years <- seq(2000, data_year - 1, 1)
        warning('data_year ', data_year, ' is not offically supported')
    } else {
        stop('data_year ', data_year, ' is not supported')
    }
    return(years)
}

#' Produce a table of forest cover change statistics for a given AOI
#'
#' For a given AOI, this function produces two tables: an annual forest loss 
#' table (in hectares, by default), and a table specifying 1) the total area of 
#' pixels that experienced forest gain and, 2) the total area of pixels that 
#' experienced both loss and gain over the full period (from 2000 through the 
#' end date of the specific product you are using, depending on the chosen 
#' \code{data_year}).  Note that forest gain and combined loss and gain are not 
#' available in the GFC product on an annualized basis.  Use 
#' \code{\link{extract_gfc}} to extract the GFC data for the AOI, and threshold 
#' it using \code{\link{threshold_gfc}} prior to running this function.
#'
#' If the\code{aoi} \code{SpatialPolygons*} object is not in the coordinate 
#' system of \code{gfc}, it will be reprojected. If there is a "label" 
#' attribute, it will be used to label the output statistics.  Otherwise, 
#' unique names ("AOI 1", "AOI 2", etc.) will be generated and used to label 
#' the output. If multiple AOIs share the same labels, statistics will be 
#' provided for the union of these AOIs.
#'
#' @seealso \code{\link{extract_gfc}}, \code{\link{threshold_gfc}}
#'
#' @export
#' @import raster
#' @import rgdal
#' @importFrom rgeos gIntersects
#' @importFrom sp spTransform CRS proj4string
#' @param aoi one or more Area of Interest (AOI) polygon(s) as a 
#' \code{SpatialPolygons*} object. See Details.
#' @param gfc extract of GFC product for a given AOI (see 
#' \code{\link{extract_gfc}}), recoded using \code{\link{threshold_gfc}}.
#' @param scale_factor how to scale the output data (from meters). Defaults to 
#' .0001 for output in hectares.
#' @param data_year which version of the Hansen data was used when
#' @return \code{list} with two elements "loss_table", a \code{data.frame} with 
#' statistics on forest loss, and "gain_table", with the area of forest gain, 
#' and area that experienced both loss and gain. The units of the output are 
#' hectares (when \code{scale_factor} is set to .0001).
gfc_stats <- function(aoi, gfc, scale_factor=.0001, data_year=2015) {
    names(gfc) <- c('forest2000', 'lossyear', 'gain', 'lossgain', 'datamask')
    gfc_boundpoly <- as(extent(gfc), 'SpatialPolygons')
    proj4string(gfc_boundpoly) <- proj4string(gfc)
    gfc_boundpoly_wgs84 <- spTransform(gfc_boundpoly, CRS('+init=epsg:4326'))
    aoi_wgs84 <- spTransform(aoi, CRS('+init=epsg:4326'))
    if (!gIntersects(gfc_boundpoly_wgs84, aoi_wgs84)) {
        stop('aoi does not intersect supplied GFC extract')
    }

    if ((((xmin(gfc) >=-180) & (xmax(gfc) <=180)) || ((xmin(gfc) >=0) & (xmax(gfc) <=360))) &&
        (ymin(gfc) >=-90) & (ymax(gfc) <= 90)) {
        # Use the included calc_pixel_area function to calculate the area of 
        # one cell in each line of the raster, allowing for accurate areal 
        # estimates of deforestation in square meters even when imagery is in 
        # WGS84.
        message('Data appears to be in latitude/longitude. Calculating cell areas on a sphere.')
        spherical_areas <- TRUE
        # Calculate the area of a single pixel in each line of the image (to 
        # avoid repeating this calculation later on)
        pixel_areas <- calc_pixel_areas(gfc)
    } else {
        spherical_areas <- FALSE
        pixel_areas <- xres(gfc) * yres(gfc)
    }

    aoi <- spTransform(aoi, CRS(proj4string(gfc)))
    if (!('label' %in% names(aoi))) {
        aoi$label <- paste('AOI', seq(1:nrow(aoi@data)))
    }

    uniq_aoi_labels <- unique(aoi$label)

    years <- gen_year_list(data_year)

    loss_table <- data.frame(year=rep(years, length(uniq_aoi_labels)),
                             aoi=rep(uniq_aoi_labels, each=length(years)))
    loss_table$cover <- 0
    loss_table$loss <- 0

    n_years <- max(years) - 2000

    gain_table <- data.frame(period=rep(paste0('2000-', max(years)), length(uniq_aoi_labels)),
                             aoi=uniq_aoi_labels,
                             gain=rep(0, length(uniq_aoi_labels)),
                             lossgain=rep(0, length(uniq_aoi_labels)))

    for (n in 1:nrow(aoi)) {
        gfc_masked <- mask(gfc, aoi[n, ], datatype='INT1U', format="GTiff", 
                           options="COMPRESS=LZW")

        # Find the first row in the loss table that applies to this aoi
        loss_table_st_row <- match(aoi[n, ]$label, loss_table$aoi)
        # No loss in first year
        loss_table$loss[loss_table_st_row] <- NA

        # Find the first row in the loss table that applies to this aoi
        gain_table_row <- match(aoi[n, ]$label, gain_table$aoi)

        # Process by blocks to avoid unnecessary reads from disk
        bs <- blockSize(gfc_masked)
        for (block_num in 1:bs$n) {
            bl_st_row <- bs$row[block_num]
            bl_nrows <- bs$nrows[block_num]
            gfc_bl <- getValuesBlock(gfc_masked, bl_st_row, bl_nrows)

            if (spherical_areas) {
                bl_pixel_areas <- rep(pixel_areas[bl_st_row:(bl_st_row + bl_nrows - 1)], each=ncol(gfc))
            } else {
                bl_pixel_areas <- pixel_areas
            }

            # Calculate initial cover, note that areas are converted to square 
            # meters using pixel size, then converted using scale_factor 
            # (default output units are hectares)
            forest2000_col <- which(names(gfc) == 'forest2000')
            loss_table$cover[loss_table_st_row] <- loss_table$cover[loss_table_st_row] +
                sum(gfc_bl[, forest2000_col] * bl_pixel_areas * scale_factor, na.rm=TRUE)

            for (i in 1:n_years) {
                lossyear_col <- which(names(gfc) == 'lossyear')
                # n + 1 because first row is year 2000, with zero loss
                loss_table$loss[loss_table_st_row + i] <- loss_table$loss[loss_table_st_row + i] +
                    sum((gfc_bl[, lossyear_col] == i) * bl_pixel_areas * scale_factor, na.rm=TRUE)
            }

            gain_col <- which(names(gfc) == 'gain')
            lossgain_col <- which(names(gfc) == 'lossgain')
            gain_table[gain_table_row, ]$gain <- gain_table[gain_table_row, ]$gain +
                sum(gfc_bl[, gain_col] * bl_pixel_areas * scale_factor, na.rm=TRUE)
            gain_table[gain_table_row, ]$lossgain <- gain_table[gain_table_row, ]$lossgain +
                sum(gfc_bl[, lossgain_col] * bl_pixel_areas * scale_factor, na.rm=TRUE)
        }

        # Calculate cover for each year by accounting for loss in prior years
        for (i in 1:n_years) {
            this_row <- loss_table_st_row + i
            loss_table$cover[this_row] <- loss_table$cover[this_row - 1] - loss_table$loss[this_row]
        }

    }

    return(list(loss_table=loss_table, gain_table=gain_table))
}
