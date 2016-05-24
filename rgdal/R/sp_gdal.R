GDALinfo <- function(fname, silent=FALSE, returnRAT=FALSE, returnCategoryNames=FALSE, returnStats=TRUE, returnColorTable=FALSE, OVERRIDE_PROJ_DATUM_WITH_TOWGS84=NULL, returnScaleOffset=TRUE) {
	if (nchar(fname) == 0) stop("empty file name")
	x <- GDAL.open(fname, silent=silent)
	d <- dim(x)[1:2]
        dr <- getDriverName(getDriver(x))
#	p4s <- .Call("RGDAL_GetProjectionRef", x, PACKAGE="rgdal")
	p4s <- getProjectionRef(x, OVERRIDE_PROJ_DATUM_WITH_TOWGS84=OVERRIDE_PROJ_DATUM_WITH_TOWGS84)
	if (nchar(p4s) == 0) p4s <- as.character(NA)
	gt <- .Call('RGDAL_GetGeoTransform', x, PACKAGE="rgdal")
        if (attr(gt, "CE_Failure") && !silent)
            warning("GeoTransform values not available")
	nbands <- .Call('RGDAL_GetRasterCount', x, PACKAGE="rgdal")
        mdata <- .Call('RGDAL_GetMetadata', x, NULL, PACKAGE="rgdal")
        subdsmdata <- .Call('RGDAL_GetMetadata', x, "SUBDATASETS",
            PACKAGE="rgdal")
        if (nbands < 1) {
#            warning("no bands in dataset")
            df <- NULL
        } else {
            band <- 1:nbands
            GDType <- character(nbands)
            hasNoDataValues <- logical(nbands)
            NoDataValues <- numeric(nbands)
            blockSize1 <- integer(nbands)
            blockSize2 <- integer(nbands)
            if (returnStats) {
              Bmin <- rep(as.numeric(NA), nbands)
              Bmax <- rep(as.numeric(NA), nbands)
              Bmn <- rep(as.numeric(NA), nbands)
              Bsd <- rep(as.numeric(NA), nbands)
            }
#            Pix <- character(nbands)
            if (returnRAT) RATlist <- vector(mode="list", length=nbands)
            if (returnCategoryNames)
                CATlist <- vector(mode="list", length=nbands)
            if (returnColorTable)
                colTabs <- vector(mode="list", length=nbands)
				
#RH 4feb2013
            if (returnScaleOffset) {
                scaleOffset <- matrix(0, ncol=2, nrow=nbands)
                colnames(scaleOffset) <- c('scale', 'offset')
            }
				
	
            for (i in seq(along = band)) {

                raster <- getRasterBand(x, band[i])
                GDType[i] <- .GDALDataTypes[(.Call("RGDAL_GetBandType",
                    raster, PACKAGE="rgdal"))+1]
                bs <- getRasterBlockSize(raster)
                blockSize1[i] <- bs[1]
                blockSize2[i] <- bs[2]
                if (returnStats) {
                  statsi <- .Call("RGDAL_GetBandStatistics", raster, silent,
                    PACKAGE="rgdal")
                  if (is.null(statsi)) {
                    Bmin[i] <- .Call("RGDAL_GetBandMinimum", raster,
                        PACKAGE="rgdal")
                    Bmax[i] <- .Call("RGDAL_GetBandMaximum", raster,
                        PACKAGE="rgdal")
                  } else {
                    Bmin[i] <- statsi[1]
                    Bmax[i] <- statsi[2]
                    Bmn[i] <- statsi[3]
                    Bsd[i] <- statsi[4]
                  }
                }
                if (returnRAT) {
                    RATi <- .Call("RGDAL_GetRAT", raster, PACKAGE="rgdal")
                    if (!is.null(RATi)) RATlist[[i]] <- RATi
                }
                if (returnCategoryNames) {
                    CATi <- .Call("RGDAL_GetCategoryNames", raster,
                        PACKAGE="rgdal")
                    if (!is.null(CATi)) CATlist[[i]] <- CATi
                }
                if (returnColorTable) {
                    colTabs[[i]] <- getBandColorTable(raster)
                }
				
#RH 4feb2013
                if (returnScaleOffset) {
                    scaleOffset[i,1] <- .Call('RGDAL_GetScale', raster,
                        PACKAGE="rgdal")
                    scaleOffset[i,2] <- .Call('RGDAL_GetOffset', raster,
                        PACKAGE="rgdal")
                }
				
                NDV <- .Call("RGDAL_GetBandNoDataValue", raster,
                    PACKAGE="rgdal")
                if (is.null(NDV)) {
                    hasNoDataValues[i] <- FALSE
                } else {
                    hasNoDataValues[i] <- TRUE
                    NoDataValues[i] <- NDV[1]
                }
#                Pix[i] <- .Call("RGDAL_GetBandMetadataItem",
#                    raster, "PIXELTYPE", "IMAGE_STRUCTURE", PACKAGE="rgdal")
            }
            df <- data.frame(GDType=GDType, hasNoDataValue=hasNoDataValues,
                NoDataValue=NoDataValues, blockSize1=blockSize1,
                blockSize2=blockSize2)
            if (returnStats) df <- cbind(df, data.frame(Bmin=Bmin,
                Bmax=Bmax, Bmean=Bmn, Bsd=Bsd))
        }
        
	GDAL.close(x)
#	res <- c(rows=d[1], columns=d[2], bands=nbands, ll.x=gt[1], ll.y=gt[4],
#		res.x=abs(gt[2]), res.y=abs(gt[6]), oblique.x=abs(gt[3]), 
#		oblique.y=abs(gt[5]))
### Modified: MDSumner 22 November 2008
        cellsize = abs(c(gt[2], gt[6]))
        ysign <- sign(gt[6])
        offset.y <- ifelse(ysign < 0, gt[4] + ysign * d[1] * abs(cellsize[2]),
            gt[4] +   abs(cellsize[2]))
        res <- c(rows = d[1], columns = d[2], bands = nbands, ll.x = gt[1],
            ll.y = offset.y, res.x = abs(gt[2]), res.y = abs(gt[6]),
            oblique.x = abs(gt[3]), oblique.y = abs(gt[5]))
#### end modification
	attr(res, "ysign") <- ysign
	attr(res, "driver") <- dr 
	attr(res, "projection") <- p4s 
	attr(res, "file") <- fname
        attr(res, "df") <- df
        attr(res, "sdf") <- returnStats
        attr(res, "mdata") <- mdata
        attr(res, "subdsmdata") <- subdsmdata
        if (returnRAT) attr(res, "RATlist") <- RATlist
        if (returnCategoryNames) attr(res, "CATlist") <- CATlist
        if (returnColorTable) attr(res, "ColorTables") <- colTabs
		
#RH 4feb2013    
        if (returnScaleOffset) attr(res, "ScaleOffset") <- scaleOffset
	class(res) <- "GDALobj"
	res
}

print.GDALobj <- function(x, ...) {
	cat("rows       ", x[1], "\n")
	cat("columns    ", x[2], "\n")
	cat("bands      ", x[3], "\n")
	cat("lower left origin.x       ", x[4], "\n")
	cat("lower left origin.y       ", x[5], "\n")
	cat("res.x      ", x[6], "\n")
	cat("res.y      ", x[7], "\n")
        cat("ysign      ", attr(x, "ysign"), "\n")
	cat("oblique.x  ", x[8], "\n")
	cat("oblique.y  ", x[9], "\n")
	cat("driver     ", attr(x, "driver"), "\n")
	cat("projection ", paste(strwrap(attr(x, "projection")),
		collapse="\n"), "\n")
	cat("file       ", attr(x, "file"), "\n")
        if (!is.null(attr(x, "df"))) {
            cat("apparent band summary:\n")
            print(attr(x, "df")[,1:5])
        }
        if (attr(x, "sdf")) {
            cat("apparent band statistics:\n")
            print(attr(x, "df")[,6:9])
        }
        if (!is.null(attr(x, "ScaleOffset"))) {
            somat <- attr(x, "ScaleOffset")
            rws <- which(somat[,1] != 1 | somat[,2] != 0)
            if (any(rws)) {
                cat("ScaleOffset:\n")
                rownames(somat) <- paste("band", 1:nrow(somat), sep="")
                print(somat[rws,])
            }
        }
        if (!is.null(attr(x, "mdata"))) {
            cat("Metadata:\n")
            cv <- attr(x, "mdata")
            for (i in 1:length(cv)) cat(cv[i], "\n")
        }
        if (!is.null(attr(x, "subdsmdata"))) {
            cat("Subdatasets:\n")
            cv <- attr(x, "subdsmdata")
            for (i in 1:length(cv)) cat(cv[i], "\n")
        }
        if (!is.null(attr(x, "RATlist"))) {
            RATs <- attr(x, "RATlist")
            nRAT <- length(RATs)
            if (nRAT == 1 ) cat("Raster attribute table:\n")
            else cat("Raster attribute tables (", nRAT, "):\n", sep="")
            for (i in 1:nRAT) {
                if (i > 1) cat("----------------------\n")
                RAT <- RATs[[i]]
                print(as.data.frame(RAT))
                cat(paste("  types:", paste(attr(RAT, "GFT_type"),
                    collapse=", ")), "\n")
                cat(paste("  usages:", paste(attr(RAT, "GFT_usage"),
                    collapse=", ")), "\n")
            }
        }
        if (!is.null(attr(x, "CATlist"))) {
            CATs <- attr(x, "CATlist")
            nCAT <- length(CATs)
            cat("Category names:\n")
            print(CATs)            
        }
        if (!is.null(attr(x, "ColorTables"))
            && length(attr(x, "ColorTables")) > 0)
            cat("Colour tables returned for bands:",
              paste(which(sapply(attr(x, "ColorTables"),
              function(x) !is.null(x))), collapse=" "), "\n")
	invisible(x)
}

asGDALROD_SGDF <- function(from) {
	x <- from
	d = dim(x)
	half.cell <- c(0.5,0.5)
	offset <- c(0,0)
	output.dim <- d[1:2]
#	p4s <- .Call("RGDAL_GetProjectionRef", x, PACKAGE="rgdal")
	p4s <- getProjectionRef(x, OVERRIDE_PROJ_DATUM_WITH_TOWGS84=NULL)
	if (nchar(p4s) == 0) p4s <- as.character(NA)
	gt = .Call('RGDAL_GetGeoTransform', x, PACKAGE="rgdal")
        if (attr(gt, "CE_Failure")) warning("GeoTransform values not available")
	if (any(gt[c(3,5)] != 0.0)) stop("Diagonal grid not permitted")
	data = getRasterData(x, list_out=TRUE)
	cellsize = abs(c(gt[2],gt[6]))
	ysign <- sign(gt[6])
	co.x <- gt[1] + (offset[2] + half.cell[2]) * cellsize[1]
	co.y <- ifelse(ysign < 0, gt[4] + (ysign*((output.dim[1] + 
		offset[1]) + (ysign*half.cell[1]))) * abs(cellsize[2]),
		gt[4] + (ysign*((offset[1]) + (ysign*half.cell[1]))) * 
		abs(cellsize[2]))
	cellcentre.offset <- c(x=co.x, y=co.y)
	grid = GridTopology(cellcentre.offset, cellsize, rev(output.dim))
#	if (length(d) == 2L)
#		df = list(band1 = as.vector(data))
#	else {
#		df <- vector(mode="list", length=d[3])
#		df[[1]] <- as.vector(data[,,1, drop = FALSE])
#		for (band in 2:d[3])
#			df[[band]] <- as.vector(data[,,band, drop = FALSE])
#		names(df) = paste("band", 1:d[3], sep="")
#	}
	return(SpatialGridDataFrame(grid = grid, 
		data = as.data.frame(data), proj4string=CRS(p4s)))
#		data = data.frame(df), proj4string=CRS(p4s)))
}

setAs("GDALReadOnlyDataset", "SpatialGridDataFrame", asGDALROD_SGDF)

asSGDF_GROD <- function(x, offset, region.dim, output.dim, p4s=NULL, ..., half.cell=c(0.5,0.5), OVERRIDE_PROJ_DATUM_WITH_TOWGS84=NULL) {
	if (!extends(class(x), "GDALReadOnlyDataset"))
		stop("x must be or extend a GDALReadOnlyDataset")
	d = dim(x)
	if (missing(offset)) offset <- c(0,0)
	if (missing(region.dim)) region.dim <- dim(x)[1:2]
	odim_flag <- NULL
	if (!missing(output.dim)) odim_flag <- TRUE
	else {
		output.dim <- region.dim
		odim_flag <- FALSE
	}

# suggestion by Paul Hiemstra 070817
	if (is.null(p4s)) 
#	    p4s <- .Call("RGDAL_GetProjectionRef", x, PACKAGE="rgdal")
	    p4s <- getProjectionRef(x, OVERRIDE_PROJ_DATUM_WITH_TOWGS84=OVERRIDE_PROJ_DATUM_WITH_TOWGS84)
	if (nchar(p4s) == 0) p4s <- as.character(NA)
	gt = .Call('RGDAL_GetGeoTransform', x, PACKAGE="rgdal")
        if (attr(gt, "CE_Failure")) warning("GeoTransform values not available")
	if (any(gt[c(3,5)] != 0.0)) stop("Diagonal grid not permitted")
	data = getRasterData(x, offset=offset, 
		region.dim=region.dim, output.dim=output.dim, ..., list_out=TRUE)
	if (!odim_flag) cellsize = abs(c(gt[2],gt[6]))
	else {
		icellsize = abs(c(gt[2],gt[6]))
		span <- icellsize * rev(d)
		cellsize <- span / rev(output.dim)
	}
	ysign <- sign(gt[6])
	co.x <- gt[1] + (offset[2] + half.cell[2]) * cellsize[1]
	co.y <- ifelse(ysign < 0, gt[4] + (ysign*((output.dim[1] + 
		offset[1]) + (ysign*half.cell[1]))) * abs(cellsize[2]),
		gt[4] + (ysign*((offset[1]) + (ysign*half.cell[1]))) * 
		abs(cellsize[2]))
	cellcentre.offset <- c(x=co.x, y=co.y)
	grid = GridTopology(cellcentre.offset, cellsize, rev(output.dim))
#	if (length(d) == 2L)
#		df = list(band1 = as.vector(data))
#	else {
#		df <- vector(mode="list", length=d[3])
#		df[[1]] <- as.vector(data[,,1, drop = FALSE])
#		for (band in 2:d[3])
#			df[[band]] <- as.vector(data[,,band, drop = FALSE])
#		names(df) = paste("band", 1:d[3], sep="")
#	}
#	df1 <- data.frame(df)
	df1 <- as.data.frame(data)
	data = SpatialGridDataFrame(grid = grid, 
		data = df1, proj4string=CRS(p4s))
	return(data)
}

readGDAL = function(fname, offset, region.dim, output.dim, band, p4s=NULL, ..., half.cell=c(0.5,0.5), silent = FALSE, OVERRIDE_PROJ_DATUM_WITH_TOWGS84=NULL) {
	if (nchar(fname) == 0) stop("empty file name")
	x = GDAL.open(fname, silent=silent)
	d = dim(x)
	if (missing(offset)) offset <- c(0,0)
	if (missing(region.dim)) region.dim <- dim(x)[1:2] # rows=nx, cols=ny
#	else d <- region.dim
	odim_flag <- NULL
        if (missing(band)) band <- NULL
        else {
                if (length(band) > 1L) d[3] <- length(band)
                else d <- d[1:2]
        }
# bug report Mike Sumner 070522
	if (!missing(output.dim)) odim_flag <- TRUE
	else {
		output.dim <- region.dim
		odim_flag <- FALSE
	}
	if (!silent) {
		cat(paste(fname, "has GDAL driver", getDriverName(getDriver(x)),"\n"))
		cat(paste("and has", d[1], "rows and", d[2], "columns\n"))
	}
# suggestion by Paul Hiemstra 070817
	if (is.null(p4s)) 
#	    p4s <- .Call("RGDAL_GetProjectionRef", x, PACKAGE="rgdal")
	    p4s <- getProjectionRef(x, OVERRIDE_PROJ_DATUM_WITH_TOWGS84=OVERRIDE_PROJ_DATUM_WITH_TOWGS84)
	if (nchar(p4s) == 0) p4s <- as.character(NA)
	gt = .Call('RGDAL_GetGeoTransform', x, PACKAGE="rgdal")
        if (attr(gt, "CE_Failure")) warning("GeoTransform values not available")
	# [1] 178400     40      0 334000      0    -40
        opSilent <- get("silent", envir=.RGDAL_CACHE)
        assign("silent", silent, envir=.RGDAL_CACHE)
	if (any(gt[c(3,5)] != 0.0)) {
		data = getRasterTable(x, band=band, offset=offset, 
			region.dim=region.dim, ...)
		GDAL.close(x)
		coordinates(data) = c(1,2)
		proj4string(data) = CRS(p4s)
	} else {
#		cellsize = abs(c(gt[2],gt[6]))
		if (!odim_flag) cellsize = abs(c(gt[2],gt[6]))
		else {
			icellsize = abs(c(gt[2],gt[6]))
# bug report Jose M. Blanco Moreno 091004
			span <- icellsize * rev(region.dim)
# bug report Mike Sumner 070215
			cellsize <- span / rev(output.dim)
		}
		ysign <- sign(gt[6])
                if (ysign > 0) 
                    warning("Y axis resolution positive, examine data for flipping")
#		cells.dim = c(d[1], d[2]) # c(d[2],d[1])
# bug report Jose M. Blanco Moreno 091004
		co.x <- gt[1] + ((offset[2]/(cellsize[1]/abs(gt[2]))) + 
                    half.cell[2]) * cellsize[1]
		co.y <- ifelse(ysign < 0, gt[4] + (ysign*((output.dim[1] + 
# bug report Jose M. Blanco Moreno 091004
			(offset[1]/(cellsize[2]/abs(gt[6]))) + 
                        (ysign*half.cell[1])))) * abs(cellsize[2]),
			gt[4] + (ysign*((offset[1]) + (ysign*half.cell[1]))) * 
			abs(cellsize[2]))
		cellcentre.offset <- c(x=co.x, y=co.y)
#		cellcentre.offset = c(x = gt[1] + 0.5 * cellsize[1], 
#			y = gt[4] - (d[2] - 0.5) * abs(cellsize[2]))
		grid = GridTopology(cellcentre.offset, cellsize, 
			rev(output.dim))
#			rev(region.dim))
		data = getRasterData(x, band=band, offset=offset, 
			region.dim=region.dim, output.dim=output.dim, ..., list_out=TRUE)
		GDAL.close(x)
#		if (length(d) == 2L)
#			df = list(band1 = as.vector(data))
#		else {
#			df <- vector(mode="list", length=d[3])
#			df[[1]] <- as.vector(data[,,1, drop = FALSE])
#			for (band in 2:d[3])
#				df[[band]] <- as.vector(data[,,band, drop = FALSE])
#			df = as.data.frame(df)
#			names(df) = paste("band", 1:d[3], sep="")
#		}
		data = SpatialGridDataFrame(grid = grid, 
			data = as.data.frame(data), proj4string=CRS(p4s))
	}
        assign("silent", opSilent, envir=.RGDAL_CACHE)
	return(data)
}

writeGDAL = function(dataset, fname, drivername = "GTiff", type = "Float32", 
		mvFlag = NA, options=NULL, copy_drivername = "GTiff",
                setStatistics=FALSE, colorTables=NULL, catNames=NULL)
{
	if (nchar(fname) == 0) stop("empty file name")
        x <- gdalDrivers()
	copy_only <- as.character(x[!x$create & x$copy, 1])
        if (drivername %in% copy_only) {
	    tds.create <- create2GDAL(dataset=dataset,
                drivername=copy_drivername, type=type,
                mvFlag=mvFlag, fname=NULL, setStatistics=setStatistics,
                colorTables=colorTables, catNames=catNames)
            tds.copy <- copyDataset(tds.create, driver=drivername, fname=fname)
            GDAL.close(tds.create)
	    saveDataset(tds.copy, fname, options=options)
# RSB 120921
            GDAL.close(tds.copy)
        } else {
	    tds.out <- create2GDAL(dataset=dataset, drivername=drivername, 
		type=type, mvFlag=mvFlag, options=options, fname=fname,
                setStatistics=setStatistics,  colorTables=colorTables,
                catNames=catNames)
	    saveDataset(tds.out, fname, options=options)
# RSB 120921
            GDAL.close(tds.out)
        }
# RSB 081030 GDAL.close cleanup
#	tmp.obj <- saveDataset(tds.out, fname, options=options)
#        GDAL.close(tmp.obj)
	invisible(fname)
}

create2GDAL = function(dataset, drivername = "GTiff", type = "Float32", mvFlag = NA, options=NULL, fname=NULL, setStatistics=FALSE, colorTables=NULL, catNames=NULL)
{
	stopifnot(gridded(dataset))
	fullgrid(dataset) = TRUE
	if (is.na(match(type, .GDALDataTypes)))
            stop(paste("Invalid type:", type, "not in:",
                paste(.GDALDataTypes, collapse="\n")))
# mvFlag issues Robert Hijmans 101109
        if (is.na(mvFlag)) {
            if (type %in% c('Byte', 'UInt16', 'Int16'))
                warning(paste("mvFlag=NA unsuitable for type", type))
        }
#	d.dim = dim(as.matrix(dataset[1])) RSB 081106
	gp = gridparameters(dataset)
	cellsize = gp$cellsize
	offset = gp$cellcentre.offset
	dims = gp$cells.dim
	d.drv = new("GDALDriver", drivername)
	nbands = length(names(slot(dataset, "data")))
        if (!is.null(options) && !is.character(options))
                stop("options not character")
	tds.out = new("GDALTransientDataset", driver = d.drv, 
		rows = dims[2], cols = dims[1],
        	bands = nbands, type = type, options = options, fname = fname,
		handle = NULL)
	gt = c(offset[1] - 0.5 * cellsize[1], cellsize[1], 0.0, 
		offset[2] + (dims[2] -0.5) * cellsize[2], 0.0, -cellsize[2])
	.Call("RGDAL_SetGeoTransform", tds.out, gt, PACKAGE = "rgdal")

	p4s <- proj4string(dataset)
	if (!is.na(p4s) && nchar(p4s) > 0) {
	    .Call("RGDAL_SetProject", tds.out, p4s, PACKAGE = "rgdal")
        } else {
            if (getDriverName(getDriver(tds.out)) == "RST") 
                stop("RST files must have a valid coordinate reference system")
        }
        if (!is.null(colorTables)) {
            stopifnot(is.list(colorTables))
            stopifnot(length(colorTables) == nbands)
            if (type != "Byte") {
#                colorTables <- NULL
                warning("colorTables valid for Byte type only in some drivers")
            }
        }
        if (!is.null(catNames)) {
            stopifnot(is.list(catNames))
            stopifnot(length(catNames) == nbands)
        }
	for (i in 1:nbands) {
		band = as.matrix(dataset[i])
		if (!is.numeric(band)) stop("Numeric bands required")
                if (setStatistics) {
                    statistics <- range(c(band), na.rm=TRUE)
                    statistics <- c(statistics, mean(c(band), na.rm=TRUE))
                    statistics <- c(statistics, sd(c(band), na.rm=TRUE))
                }
		if (!is.na(mvFlag))
			band[is.na(band)] = mvFlag
		putRasterData(tds.out, band, i)
                tds.out_b <- getRasterBand(dataset=tds.out, band=i)
		if (!is.na(mvFlag)) {
		    .Call("RGDAL_SetNoDataValue", tds.out_b, as.double(mvFlag),
		        PACKAGE = "rgdal")
		}
                if (setStatistics) {
                    .gd_SetStatistics(tds.out_b, as.double(statistics))
                }
                if (!is.null(colorTables)) {
                    icT <- colorTables[[i]]
                    if (!is.null(icT)) {
                        .gd_SetRasterColorTable(tds.out_b, icT)
                    }
                }
                if (!is.null(catNames)) {
                    icN <- catNames[[i]]
                    if (!is.null(icN)) {
                        .gd_SetCategoryNames(tds.out_b, icN)
                    }       
                }
	}
	tds.out
}


gdalDrivers <- function() getGDALDriverNames()

toSigned <- function(x, base) {
    if (any(x < 0)) stop("already signed")
    if (storage.mode(x) != "integer") stop("band not integer")
    b_2 <- (2^(base-1)-1)
    b <- 2^base
    x[x > b_2] <- x[x > b_2] - b
    as.integer(x)
}

toUnSigned <- function(x, base) {
    if (all(x >= 0)) stop("already unsigned")
    if (storage.mode(x) != "integer") stop("band not integer")
    b <- 2^base
    x[x < 0] <- x[x < 0] + b
    as.integer(x)
}

"GDALSpatialRef" <- function(fname, silent=FALSE, OVERRIDE_PROJ_DATUM_WITH_TOWGS84=NULL) {
	if (nchar(fname) == 0) stop("empty file name")
	x <- GDAL.open(fname, silent=silent)
#        p4s <- .Call("RGDAL_GetProjectionRef", x, PACKAGE="rgdal")
        p4s <- getProjectionRef(x, OVERRIDE_PROJ_DATUM_WITH_TOWGS84=OVERRIDE_PROJ_DATUM_WITH_TOWGS84)
	GDAL.close(x)
        p4s
}

