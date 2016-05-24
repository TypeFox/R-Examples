SGDF2PCT <- function(x, ncolors=256, adjust.bands=TRUE) {
	if (!inherits(x, "SpatialGridDataFrame")) 
		stop("SpatialGridDataFrame required")
	if (length(names(slot(x, "data"))) != 3L)
		stop("Three data columns (red, green, blue) required")
	if (ncolors > 256) {
		warning("ncolors reset to maximum 256")
		ncolors <- 256
	}
	if (ncolors < 2) {
		warning("ncolors reset to minimum 2")
		ncolors <- 2
	}

	fullgrid(x) = TRUE
	d.dim <- dim(as.matrix(x[1]))
	d.drv <- new("GDALDriver", "GTiff")
	GTiff3B <- new("GDALTransientDataset", driver = d.drv, 
		rows = d.dim[2], cols = d.dim[1], bands = 3, type = "Byte", 
		handle = NULL)
	gp <- gridparameters(x)
	cellsize <- gp$cellsize
	offset <- gp$cellcentre.offset
	dims <- gp$cells.dim
	gt <- c(offset[1] - 0.5 * cellsize[1], cellsize[1], 0.0, 
		offset[2] + (dims[2] -0.5) * cellsize[2], 0.0, -cellsize[2])
	.Call("RGDAL_SetGeoTransform", GTiff3B, gt, PACKAGE = "rgdal")
	have_NAs <- FALSE
	for (i in 1:3) {
		band = as.matrix(x[i])
		if (any(is.na(band))) have_NAs <- TRUE
		if (!is.numeric(band)) stop("Numeric bands required")
# 101213 Michael Sumner
		if (adjust.bands || have_NAs) {
			bmax <- max(band, na.rm = TRUE)
			bmin <- min(band, na.rm = TRUE)
			if (bmax == bmin) {
                            if (ncolors < 256) bmax <- bmin + ncolors
                            else bmax <- bmin + 1
                        }
			band <- floor((band - bmin)/((bmax-bmin)/(255)))
		} else {
			if (!is.integer(band)) 
				stop("unadjusted band not integer")
			if (any(band[!is.na(band)] < 0) || 
				any(band[!is.na(band)] > 255)) 
				stop("unadjusted band out of range")
		}
		storage.mode(band) <- "integer"
		putRasterData(GTiff3B, band, i)
	}
# 101213 Michael Sumner
#	if (have_NAs) ncolors <- ncolors + 1
	dx <- RGB2PCT(GTiff3B, band=1:3, ncolors=ncolors, set.ctab=FALSE)
	GDAL.close(GTiff3B)
	output <- getRasterData(dx$dataset, band=1)
	idx <- as.vector(output)+1
	ct <- dx$pct
	res <- list(idx=idx, ct=ct[1:ncolors])
	GDAL.close(dx$dataset)
	res
}

vec2RGB <- function(vec, breaks, col) {
	if (!is.numeric(vec)) stop("vec must be numeric")
	if (!is.vector(vec)) stop("vec must be a vector")
	if (length(col) != (length(breaks)-1L)) 
		stop("length of col must be one less than length of breaks")
	idvec <- findInterval(vec, breaks, all.inside=TRUE)
	rgb_col <- col2rgb(col)
	res <- rgb_col[, idvec]
	storage.mode(res) <- "integer"
	t(res)
}
