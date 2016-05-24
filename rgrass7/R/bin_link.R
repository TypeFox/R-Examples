# Interpreted GRASS 7 interface functions
# Copyright (c) 2015 Roger S. Bivand
#

readRAST <- function(vname, cat=NULL, ignore.stderr = get.ignore.stderrOption(), 
	NODATA=NULL, plugin=get.pluginOption(), mapset=NULL, useGDAL=get.useGDALOption(), close_OK=TRUE,
        drivername="GTiff", driverFileExt=NULL, return_SGDF=TRUE) {
	if (!is.null(cat))
		if(length(vname) != length(cat)) 
			stop("vname and cat not same length")
    if (get.suppressEchoCmdInFuncOption()) {
        inEchoCmd <- set.echoCmdOption(FALSE)
    }
                if (close_OK) {
                    openedConns <- as.integer(row.names(showConnections()))
                }

        tryCatch(
            {
                stopifnot(is.logical(plugin)|| is.null(plugin))
                if (!is.null(plugin) && plugin && length(vname) > 1) plugin <- FALSE
                stopifnot(is.logical(ignore.stderr))
                stopifnot(is.logical(useGDAL))
                if (useGDAL) {
                    if (requireNamespace("rgdal", quietly = TRUE)) {
                        gdalD <- rgdal::gdalDrivers()$name
                    } else {
                        stop("rgdal not available")
                    }
                }
                if (!useGDAL && is.null(plugin)) plugin <- FALSE
                
                if (is.null(plugin)) plugin <- "GRASS" %in% gdalD
                if (length(vname) > 1) plugin <- FALSE
                if (plugin) {
                    resa <- .read_rast_plugin(vname, mapset=mapset, ignore.stderr=ignore.stderr)
                } else {
                    resa <- .read_rast_non_plugin(vname=vname, NODATA=NODATA,
                                                  driverFileExt=driverFileExt, ignore.stderr=ignore.stderr, return_SGDF=return_SGDF, cat=cat)
                }
            },
            finally = {
                if (close_OK) { #closeAllConnections()
                    openConns_now <- as.integer(row.names(showConnections()))
                    toBeClosed <- openConns_now[!(openConns_now %in% openedConns)]
                    for (bye in toBeClosed) close(bye)
                }
                if (get.suppressEchoCmdInFuncOption()) {
                    tull <- set.echoCmdOption(inEchoCmd)
                }
            }
        )

    resa
}


.read_rast_non_plugin <- function(vname, NODATA, driverFileExt, ignore.stderr, return_SGDF, cat){
    {
	pid <- as.integer(round(runif(1, 1, 1000)))
	p4 <- CRS(getLocationProj())

        reslist <- vector(mode="list", length=length(vname))
        names(reslist) <- vname

	for (i in seq(along=vname)) {


                glist <- execGRASS("r.info", flags="g", map=vname[i],
                               intern=TRUE, ignore.stderr=ignore.stderr)
                whCELL <- glist[grep("datatype", glist)]
		to_int <- length(which(unlist(strsplit(
			whCELL, "=")) == "CELL")) > 0
                Dcell <- length(which(unlist(strsplit(
			whCELL, "=")) == "DCELL")) > 0

		gtmpfl1 <- dirname(execGRASS("g.tempfile",
		    pid=pid, intern=TRUE,
		    ignore.stderr=ignore.stderr))
		rtmpfl1 <- ifelse(.Platform$OS.type == "windows" &&
                        (Sys.getenv("OSTYPE") == "cygwin"), 
			system(paste("cygpath -w", gtmpfl1, sep=" "), 
			intern=TRUE), gtmpfl1)
		gtmpfl11 <- paste(gtmpfl1, vname[i], sep=.Platform$file.sep)
		rtmpfl11 <- paste(rtmpfl1, vname[i], sep=.Platform$file.sep)
                if (!is.null(driverFileExt)) {
                    gtmpfl11 <- paste(gtmpfl1, driverFileExt, sep=".")
                    rtmpfl11 <- paste(rtmpfl1, driverFileExt, sep=".")
                }

                if (!is.null(NODATA)) {
		    if (!is.finite(NODATA) || !is.numeric(NODATA))
			stop("invalid NODATA value")
		    if (NODATA != round(NODATA)) 
			warning("NODATA rounded to integer")
		    NODATA <- round(NODATA)
		}

# 130422 at rgdal 0.8-8 GDAL.close(DS)
# 061107 Dylan Beaudette NODATA
# 071009 Markus Neteler's idea to use range
	  	    if (is.null(NODATA)) {
		      tx <- execGRASS("r.info", flags="r",
		        map=vname[i], intern=TRUE, 
		        ignore.stderr=ignore.stderr)
		      tx <- gsub("=", ":", tx)
		      con <- textConnection(tx)
		      res <- read.dcf(con)
		      close(con)
		      lres <- as.list(res)
		      names(lres) <- colnames(res)
		      lres <- lapply(lres, function(x)
                          ifelse(x == "NULL", as.numeric(NA), as.numeric(x)))
		      if (!is.numeric(lres$min) || 
			!is.finite(as.double(lres$min))) 
			    NODATA <- as.integer(999)
		      else {
			lres$min <- floor(as.double(lres$min))
		        NODATA <- floor(lres$min) - 1
		      }
                    }
                    rOutBinFlags <- "b"
# 120118 Rainer Krug
                        if (to_int) rOutBinFlags <- c(rOutBinFlags, "i")
                        else rOutBinFlags <- c(rOutBinFlags, "f")
                        rOutBinBytes <- 4L
                if (Dcell) rOutBinBytes <- 8L
                tryCatch(
                    {
		        execGRASS("r.out.bin", flags=rOutBinFlags, 
                                  input=vname[i], output=gtmpfl11, bytes=rOutBinBytes,
                                  null=as.integer(NODATA),
                                  ignore.stderr=ignore.stderr)
                        
                        gdal_info <- bin_gdal_info(rtmpfl11, to_int)
                        
                        what <- ifelse(to_int, "integer", "double")
                        n <- gdal_info[1] * gdal_info[2]
                        size <- gdal_info[10]/8
                        
                        reslist[[i]] <- readBinGridData(rtmpfl11, what=what,
                                                        n=n, size=size, endian=attr(gdal_info, "endian"),
                                                        nodata=attr(gdal_info, "nodata"))
                    },
                    finally = {
                        unlink(paste(rtmpfl1, list.files(rtmpfl1,
                                                         pattern=vname[i]), sep=.Platform$file.sep))
                    }
                )

#		if (i == 1) resa <- res
#		else {
#			grida <- getGridTopology(resa)
#			grid <- getGridTopology(res)
#			if (!isTRUE(all.equal(grida, grid)))
#				stop("topology is not equal")
#			onames <- c(names(resa@data), names(res@data))
#			ncols <- dim(resa@data)[2]
#			lst <- vector(mode="list", length=ncols+1)
#			names(lst) <- onames
#			for (i in 1:ncols) lst[[i]] <- resa@data[[i]]
#			lst[[ncols+1]] <- res@data[[1]]
#			df <- data.frame(lst)
#			resa <- SpatialGridDataFrame(grid=grida, 
#				data=df, proj4string=p4)
#		}

	}


        co <- unname(c((gdal_info[4] + (gdal_info[6]/2)),
            (gdal_info[5] + (gdal_info[7]/2))))
	grid <- GridTopology(co, unname(c(gdal_info[6], gdal_info[7])),
            unname(c(gdal_info[2], gdal_info[1])))

        if (!return_SGDF) {
           res <- list(grid=grid, dataList=reslist, proj4string=p4)
           class(res) <- "gridList"
           return(res)
        }

        if (length(unique(sapply(reslist, length))) != 1)
            stop ("bands differ in length")

        df <- as.data.frame(reslist)

	resa <- SpatialGridDataFrame(grid=grid, data=df, proj4string=p4)

	if (!is.null(cat)) {
		for (i in seq(along=cat)) {
			if (cat[i] && is.integer(resa@data[[i]])) {

				rSTATS <- execGRASS("r.stats",
				    flags=c("l", "quiet"),
				    input=vname[i],
				    intern=TRUE, ignore.stderr=ignore.stderr)

				cats <- strsplit(rSTATS, " ")
				catnos <- sapply(cats, function(x) x[1])
				catlabs <- sapply(cats, 
					function(x) paste(x[-1], collapse=" "))
				if (any(!is.na(match(catnos, "*")))) {
					isNA <- which(catnos == "*")
					catnos <- catnos[-isNA]
					catlabs <- catlabs[-isNA]
				}
                                if (length(catlabs) > length(unique(catlabs))) {
                                    catlabs <- paste(catlabs, catnos, sep="_")
                                    warning("non-unique category labels; category number appended")
                                }
				resa@data[[i]] <- factor(resa@data[[i]], 
					levels=catnos, labels=catlabs)
			}
		}
	} 
    }
    return(resa)
}




bin_gdal_info <- function(fname, to_int) {
	if (!file.exists(fname)) stop(paste("no such file:", fname))
	if (!file.exists(paste(fname, "hdr", sep="."))) 
		stop(paste("no such file:", paste(fname, "hdr", sep=".")))
	if (!file.exists(paste(fname, "wld", sep="."))) 
		stop(paste("no such file:", paste(fname, "wld", sep=".")))
	con <- file(paste(fname, "hdr", sep="."), "r")
	l8 <- readLines(con, n=8)
	close(con)
	l8 <- read.dcf(con <- textConnection(gsub(" ", ":", l8)))
	close(con)
	lres <- as.list(l8)
	names(lres) <- colnames(l8)
	lres$nrows <- as.integer(lres$nrows)
	lres$ncols <- as.integer(lres$ncols)
	lres$nbands <- as.integer(lres$nbands)
	lres$nbits <- as.integer(lres$nbits)
	lres$skipbytes <- as.integer(lres$skipbytes)
	lres$nodata <- ifelse(to_int, as.integer(lres$nodata), 
		as.numeric(lres$nodata))
	lres$byteorder <- as.character(lres$byteorder)
	endian <- .Platform$endian
	if ((endian == "little" && lres$byteorder == "M") ||
		(endian == "big" && lres$byteorder == "I")) endian <- "swap"
	con <- file(paste(fname, "wld", sep="."), "r")
	l6 <- readLines(con, n=6)
	close(con)
	lres$ewres <- abs(as.numeric(l6[1]))
	lres$nsres <- abs(as.numeric(l6[4]))
	lres$n_cc <- as.numeric(l6[6])
	lres$w_cc <- as.numeric(l6[5])
	lres$s_cc <- lres$n_cc - lres$nsres * (lres$nrows-1)

    outres <- numeric(10)
    outres[1] <- lres$nrows
    outres[2] <- lres$ncols
    outres[3] <- lres$nbands
    outres[4] <- lres$w_cc - (lres$ewres/2)
    outres[5] <- lres$s_cc - (lres$nsres/2)
    outres[6] <- lres$ewres
    outres[7] <- lres$nsres
    outres[10] <- lres$nbits
    attr(outres, "endian") <- endian
    attr(outres, "nodata") <- lres$nodata
#	grid = GridTopology(c(lres$w_cc, lres$s_cc), 
#		c(lres$ewres, lres$nsres), c(lres$ncols,lres$nrows))
    outres
}

readBinGridData <- function(fname, what, n, size, endian, nodata) {
	map <- readBin(fname, what=what, n=n, size=size, signed=TRUE,
		endian=endian)
	is.na(map) <- map == nodata
	map
}

#readBinGrid <- function(fname, colname=basename(fname), 
#	proj4string=CRS(as.character(NA)), integer) {
#	if (missing(integer)) stop("integer TRUE/FALSE required")
#	if (!file.exists(fname)) stop(paste("no such file:", fname))
#	if (!file.exists(paste(fname, "hdr", sep="."))) 
#		stop(paste("no such file:", paste(fname, "hdr", sep=".")))
#	if (!file.exists(paste(fname, "wld", sep="."))) 
#		stop(paste("no such file:", paste(fname, "wld", sep=".")))
#	con <- file(paste(fname, "hdr", sep="."), "r")
#	l8 <- readLines(con, n=8)
#	close(con)
#	l8 <- read.dcf(con <- textConnection(gsub(" ", ":", l8)))
#	close(con)
#	lres <- as.list(l8)
#	names(lres) <- colnames(l8)
#	lres$nrows <- as.integer(lres$nrows)
#	lres$ncols <- as.integer(lres$ncols)
#	lres$nbands <- as.integer(lres$nbands)
#	lres$nbits <- as.integer(lres$nbits)
#	lres$skipbytes <- as.integer(lres$skipbytes)
#	lres$nodata <- ifelse(integer, as.integer(lres$nodata), 
#		as.numeric(lres$nodata))
#	lres$byteorder <- as.character(lres$byteorder)
#	endian <- .Platform$endian
#	if ((endian == "little" && lres$byteorder == "M") ||
#		(endian == "big" && lres$byteorder == "I")) endian <- "swap"
#	con <- file(paste(fname, "wld", sep="."), "r")
#	l6 <- readLines(con, n=6)
#	close(con)
#	lres$ewres <- abs(as.numeric(l6[1]))
#	lres$nsres <- abs(as.numeric(l6[4]))
#	lres$n_cc <- as.numeric(l6[6])
#	lres$w_cc <- as.numeric(l6[5])
#	lres$s_cc <- lres$n_cc - lres$nsres * (lres$nrows-1)
#
#	what <- ifelse(integer, "integer", "double")
#	n <- lres$nrows * lres$ncols
#	size <- lres$nbits/8
#	map <- readBin(fname, what=what, n=n, size=size, signed=TRUE,
#		endian=endian)
#	is.na(map) <- map == lres$nodata
#	grid = GridTopology(c(lres$w_cc, lres$s_cc), 
#		c(lres$ewres, lres$nsres), c(lres$ncols,lres$nrows))
#	df <- list(var1=map)
#	names(df) <- colname
#	df1 <- data.frame(df)
#
#	res <- SpatialGridDataFrame(grid, data = df1, proj4string=proj4string)
#	res
#}

.read_rast_plugin <- function(vname, mapset=NULL, ignore.stderr=NULL) {

        if (is.null(ignore.stderr))
            ignore.stderr <- get.ignore.stderrOption()
        stopifnot(is.logical(ignore.stderr))
        if (length(vname) > 1) stop("single raster required for plugin")

        gg <- gmeta()
        if (is.null(mapset)) {
            c_at <- strsplit(vname[1], "@")[[1]]
            if (length(c_at) == 1) {
                mapset <- .g_findfile(vname[1], type="cell")
            } else if (length(c_at) == 2) {
                mapset <- c_at[2]
                vname[1] <- c_at[1]
            } else stop("malformed raster name")
        }

        fname <- paste(gg$GISDBASE, gg$LOCATION_NAME, mapset,
            "cellhd", vname[1], sep="/")
      if (requireNamespace("rgdal", quietly = TRUE)) {
        fninfo <- rgdal::GDALinfo(fname, silent=ignore.stderr)
        chks <- logical(4)
        names(chks) <- c("cols", "rows", "origin.northing",
            "origin.easting")
        chks[1] <- isTRUE(all.equal(abs((gg$w-gg$e)/gg$ewres), fninfo[2],
            tol=2e-7, check.attributes=FALSE))
        chks[2] <- isTRUE(all.equal(abs((gg$n-gg$s)/gg$nsres), fninfo[1],
            tol=2e-7, check.attributes=FALSE))
# changed from gg$n 100129, thanks to Rainer Krug
        chks[3] <- isTRUE(all.equal(gg$s, fninfo[5],
            check.attributes=FALSE))
        chks[4] <- isTRUE(all.equal(gg$w, fninfo[4],
            check.attributes=FALSE))
        if (any(!chks)) {
            cat("raster map/current region mismatch detected in components:\n")
            print(chks)
	    stop("rerun with plugin=FALSE\n") 
        }

        resa <- rgdal::readGDAL(fname, silent=ignore.stderr)
	names(resa) <- make.names(vname)
        resa
      } else {
        stop("rgdal not available")
      }
}
    

writeRAST <- function(x, vname, zcol = 1, NODATA=NULL, 
	ignore.stderr = get.ignore.stderrOption(), useGDAL=get.useGDALOption(), overwrite=FALSE, flags=NULL,
        drivername="GTiff") {

        if (get.suppressEchoCmdInFuncOption()) {
            inEchoCmd <- set.echoCmdOption(FALSE)
        }

        tryCatch(
            {
                stopifnot(is.logical(ignore.stderr))
                stopifnot(is.logical(useGDAL))
                pid <- as.integer(round(runif(1, 1, 1000)))
                gtmpfl1 <- dirname(execGRASS("g.tempfile", pid=pid,
                                             intern=TRUE, ignore.stderr=ignore.stderr))
                
                rtmpfl1 <- ifelse(.Platform$OS.type == "windows" &&
                                      (Sys.getenv("OSTYPE") == "cygwin"), 
                                  system(paste("cygpath -w", gtmpfl1, sep=" "), intern=TRUE), 
                                  gtmpfl1)
                
                fid <- paste("X", pid, sep="")
                gtmpfl11 <- paste(gtmpfl1, fid, sep=.Platform$file.sep)
                rtmpfl11 <- paste(rtmpfl1, fid, sep=.Platform$file.sep)
                if (!is.numeric(x@data[[zcol]])) 
                    stop("only numeric columns may be exported")
                if (overwrite && !("overwrite" %in% flags))
                    flags <- c(flags, "overwrite")
                tryCatch(
                    {
                        res <- writeBinGrid(x, rtmpfl11, attr = zcol, na.value = NODATA)
                        
                        flags <- c(res$flag, flags)
                        
                        execGRASS("r.in.bin", flags=flags,
                                  input=gtmpfl11,
                                  output=vname, bytes=as.integer(res$bytes), 
                                  north=as.numeric(res$north), south=as.numeric(res$south), 
                                  east=as.numeric(res$east), west=as.numeric(res$west), 
                                  rows=as.integer(res$rows), cols=as.integer(res$cols), 
                                  anull=as.numeric(res$anull), ignore.stderr=ignore.stderr)
                    },
                    finally = {
                        unlink(paste(rtmpfl1, list.files(rtmpfl1, pattern=fid), 
                                     sep=.Platform$file.sep))
                    }
                )
            },
            finally = {
                if (get.suppressEchoCmdInFuncOption()) {
                    tull <- set.echoCmdOption(inEchoCmd)
                }
            }
        )

	invisible(res)
}

writeBinGrid <- function(x, fname, attr = 1, na.value = NULL) { 
	if (!gridded(x))
		stop("can only write SpatialGridDataFrame objects to binary grid")
	x = as(x, "SpatialGridDataFrame")
	gp = gridparameters(x)
	if (length(gp$cells.dim) != 2)
		stop("binary grid only supports 2D grids")
	z = x@data[[attr]]
	if (is.factor(z)) z <- as.numeric(z)
	if (!is.numeric(z)) stop("only numeric values may be exported")
	if (is.null(na.value)) {
		na.value <- floor(min(z, na.rm=TRUE)) - 1
	} else {
		if (!is.finite(na.value) || !is.numeric(na.value))
			stop("invalid NODATA value")
		if (na.value != round(na.value)) 
			warning("NODATA rounded to integer")
		na.value <- round(na.value)
	}
	res <- list()
	res$anull <- formatC(na.value, format="d")
	z[is.na(z)] = as.integer(na.value)
	if (storage.mode(z) == "integer") {
		sz <- 4
	} else if (storage.mode(z) == "double") {
		sz <- 8
		res$flag <- "d"
	} else stop("unknown storage mode")
	res$bytes <- formatC(sz, format="d")
	f = file(fname, open = "wb")
	writeBin(z, con=f, size=sz)
	close(f)
	grd <- slot(x, "grid")
	f = file(paste(fname, "hdr", sep="."), open = "wt")
	writeLines(paste("nrows", grd@cells.dim[2]), f)
	res$rows <- formatC(grd@cells.dim[2], format="d")
	writeLines(paste("ncols", grd@cells.dim[1]), f)
	res$cols <- formatC(grd@cells.dim[1], format="d")
	writeLines(paste("nbands 1"), f)
	writeLines(paste("nbits", 8*sz), f)
	writeLines(paste("byteorder", ifelse(.Platform$endian == "little", 
		"I", "M")), f)
	writeLines(paste("layout bil"), f)
	writeLines(paste("skipbytes 0"), f)
	writeLines(paste("nodata", na.value), f)
	close(f)
	f = file(paste(fname, "wld", sep="."), open = "wt")
	writeLines(formatC(grd@cellsize[1], format="f"), f)
	writeLines("0.0", f)
	writeLines("0.0", f)
	writeLines(formatC(-grd@cellsize[2], format="f"), f)
	writeLines(formatC(grd@cellcentre.offset[1], format="f"), f)
	writeLines(formatC((grd@cellcentre.offset[2] +
		grd@cellsize[2]*(grd@cells.dim[2]-1)), format="f"), f)
	close(f)
	res$north <- formatC((grd@cellcentre.offset[2] +
		grd@cellsize[2]*(grd@cells.dim[2]-1)
		+ 0.5*grd@cellsize[2]), format="f")
	res$south <- formatC(grd@cellcentre.offset[2] - 0.5*grd@cellsize[2], 
		format="f")
	res$east <- formatC((grd@cellcentre.offset[1] + 
		grd@cellsize[1]*(grd@cells.dim[1]-1) + 0.5*grd@cellsize[1]), 
		format="f")
	res$west <- formatC(grd@cellcentre.offset[1] - 0.5*grd@cellsize[1], 
		format="f")
	invisible(res)
}

