# Interpreted GRASS 6+ interface functions
# Copyright (c) 2005-2012 Roger S. Bivand
#

readVECT <- readVECT6 <- function(vname, layer, type=NULL, plugin=NULL,
        remove.duplicates=TRUE, 
	ignore.stderr = NULL, with_prj=TRUE, with_c=FALSE, mapset=NULL, 
	pointDropZ=FALSE, driver="ESRI Shapefile") {

    if (get.suppressEchoCmdInFuncOption()) {
        inEchoCmd <- get.echoCmdOption()
        tull <- set.echoCmdOption(FALSE)
    }
    if (is.null(plugin))
        plugin <- get("plugin", envir = .GRASS_CACHE)
    stopifnot(is.logical(plugin)|| is.null(plugin))
    if (is.null(ignore.stderr))
        ignore.stderr <- get("ignore.stderr", envir = .GRASS_CACHE)
    stopifnot(is.logical(ignore.stderr))
    if (missing(layer)) layer <- 1L
# 120908 emails Markus Neteler, Markus Metz, default TRUE before G7
    stopifnot(is.logical(with_c))
    with_c <- !with_c
    if (!ignore.stderr) 
      message("with_c: argument reversed from version 0.7-11 and in GRASS 6")

    if (driver == "GRASS") plugin <- TRUE

    if (requireNamespace("rgdal", quietly = TRUE)) {
    if (is.null(plugin)) {
        ogrD <- rgdal::ogrDrivers()$name
	plugin <- "GRASS" %in% ogrD
    }
    sss <- strsplit(packageDescription("rgdal")$Version, "-")[[1]]
    if (plugin) {
        ogrD <- rgdal::ogrDrivers()$name
	if (!("GRASS" %in% ogrD)) stop("no GRASS plugin driver")
        gg <- gmeta()
        if (is.null(mapset)) {
            c_at <- strsplit(vname[1], "@")[[1]]
            if (length(c_at) == 1) {
                mapset <- .g_findfile(vname[1], type="vector")
            } else if (length(c_at) == 2) {
                mapset <- c_at[2]
                vname[1] <- c_at[1]
            } else stop("malformed vector name")
        }
        dsn <- paste(gg$GISDBASE, gg$LOCATION_NAME, mapset,
            "vector", vname[1], "head", sep="/")
        if (sss[1] >= "0.6" && as.integer(sss[2]) > 7) {
	    res <- rgdal::readOGR(dsn, layer=as.character(layer),
                verbose=!ignore.stderr, pointDropZ=pointDropZ)
        } else {
	    res <- rgdal::readOGR(dsn, layer=as.character(layer),
                verbose=!ignore.stderr)
        }
    } else {
        ogrD <- rgdal::ogrDrivers()$name
	if (!(driver %in% ogrD))
            stop(paste("Requested driver", driver, "not available in rgdal"))
        ogrDGRASS <- execGRASS("v.in.ogr", flags="f", intern=TRUE,
            ignore.stderr=ignore.stderr)
        ogrDGRASSs <- strsplit(ogrDGRASS, ": ")
        if (!(driver %in% sapply(ogrDGRASSs, "[", 2)))
            stop(paste("Requested driver", driver, "not available in GRASS"))
        fDrivers <- c("GML", "SQLite")
        dDrivers <- c("ESRI_Shapefile", "MapInfo_File")
        if (!(gsub(" ", "_", driver) %in% c(fDrivers, dDrivers)))
            stop(paste("Requested driver", driver, "not supported"))
        is_dDriver <- TRUE
        if (gsub(" ", "_", driver) %in% fDrivers) is_dDriver <- FALSE
	vinfo <- vInfo(vname)
	types <- names(vinfo)[which(vinfo > 0)]
	if (is.null(type)) {
		if (length(grep("points", types)) > 0) type <- "point"
		if (length(grep("lines", types)) > 0) type <- "line"
		if (length(grep("areas", types)) > 0) type <- "area"
		if (is.null(type)) stop("Vector type not found")
	}

	pid <- as.integer(round(runif(1, 1, 1000)))

	gtmpfl1 <- dirname(execGRASS("g.tempfile", pid=pid,
            intern=TRUE, ignore.stderr=ignore.stderr))
	rtmpfl1 <- ifelse(.Platform$OS.type == "windows" &&
                (Sys.getenv("OSTYPE") == "cygwin"), 
		system(paste("cygpath -w", gtmpfl1, sep=" "), intern=TRUE), 
		gtmpfl1)

	if (driver == "ESRI Shapefile") 
	        shname <- substring(vname, 1, ifelse(nchar(vname) > 8, 8, 
		nchar(vname)))

        flags <- NULL
        if (with_prj) flags <- "e"
        if (with_c) flags <- c(flags, "c")
        if (is_dDriver){
            GDSN <- gtmpfl1
            RDSN <- rtmpfl1
            LAYER <- shname
        } else {
            GDSN <- paste(gtmpfl1, shname, sep=.Platform$file.sep)
            RDSN <- paste(rtmpfl1, shname, sep=.Platform$file.sep)
            LAYER <- shname
        }
        tryCatch(
            {
                execGRASS("v.out.ogr", flags=flags, input=vname,
                          type=type, layer=layer, dsn=GDSN, olayer=LAYER,
                          format=gsub(" ", "_", driver), ignore.stderr=ignore.stderr)
                
                if (sss[1] >= "0.6" && as.integer(sss[2]) > 7) {
                    res <- rgdal::readOGR(dsn=RDSN, layer=LAYER, verbose=!ignore.stderr, 
                                          pointDropZ=pointDropZ)
                } else {
                    res <- rgdal::readOGR(dsn=rtmpfl1, layer=shname, verbose=!ignore.stderr)           }
            },
            finally = {
                if (.Platform$OS.type != "windows") {
                    unlink(paste(rtmpfl1, list.files(rtmpfl1, pattern=shname), 
                                 sep=.Platform$file.sep))
                }
            }
        )
        
	if (remove.duplicates && type != "point") {
		dups <- duplicated(slot(res, "data"))
		if (any(dups)) {
			if (length(grep("line", type)) > 0) type <- "line"
			if (length(grep("area", type)) > 0) type <- "area"
			if (type != "area" && type != "line")
			    stop("try remove.duplicates=FALSE")
			ndata <- as(res, "data.frame")[!dups,,drop=FALSE]
                        cand <- as.character(ndata$cat)
                        cand[is.na(cand)] <- "na"
			row.names(ndata) <- cand
			if (type == "area") {
				pls <- slot(res, "polygons")
			} else if (type == "line") {
				pls <- slot(res, "lines")
			}
			p4s <- proj4string(res)
			IDs <- as.character(res$cat)
			IDs[is.na(IDs)] <- "na"
			tab <- table(factor(IDs))
			n <- length(tab)
			if (n + sum(dups) != length(pls))
				stop("length mismatch in duplicate removal")
			IDss <- .mixedsort(names(tab))
			reg <- match(IDs, IDss)
			belongs <- lapply(1:n, function(x) which(x == reg))
			npls <- vector(mode="list", length=n)
			for (i in 1:n) {
				nParts <- length(belongs[[i]])
				srl <- NULL
				for (j in 1:nParts) {
					plij <- pls[[belongs[[i]][j]]]
					if (type == "area") {
						plijp <- slot(plij, "Polygons")
					} else if (type == "line") {
						plijp <- slot(plij, "Lines")
					}
					srl <- c(srl, plijp)
				}
				if (type == "area") {
				    npls[[i]] <- Polygons(srl, ID=IDss[i])
				} else if (type == "line") {
				    npls[[i]] <- Lines(srl, ID=IDss[i])
				}
			}
			if (type == "area") {
			    SP <- SpatialPolygons(npls, proj4string=CRS(p4s))
			    res <- SpatialPolygonsDataFrame(SP, ndata)
			} else if (type == "line") {
			    SP <- SpatialLines(npls, proj4string=CRS(p4s))
			    res <- SpatialLinesDataFrame(SP, ndata)
			}
		}

	}
    }
    if (get.suppressEchoCmdInFuncOption()) {
        tull <- set.echoCmdOption(inEchoCmd)
    }
    } else {
      stop("rgdal not available")
    }

    res
}

# Function mixedorder copied from gtools 2.2.3 LGPL Gregory R. Warnes
.mixedsort <- function (x) {
    x[.mixedorder(x)]
}

.mixedorder <- function (x) {
    delim = "\\$\\@\\$"
    numeric <- function(x) {
        optwarn = options("warn")
        on.exit(options(optwarn))
        options(warn = -1)
        as.numeric(x)
    }
    nonnumeric <- function(x) {
        optwarn = options("warn")
        on.exit(options(optwarn))
        options(warn = -1)
        ifelse(is.na(as.numeric(x)), toupper(x), NA)
    }
    x <- as.character(x)
    which.nas <- which(is.na(x))
    which.blanks <- which(x == "")
    if (length(which.blanks) > 0) 
        x[which.blanks] <- -Inf
    if (length(which.nas) > 0) 
        x[which.nas] <- Inf
    delimited <- gsub("([+-]{0,1}[0-9.]+([eE][+-]{0,1}[0-9.]+){0,1})", 
        paste(delim, "\\1", delim, sep = ""), x)
    step1 <- strsplit(delimited, delim)
    step1 <- lapply(step1, function(x) x[x > ""])
    step1.numeric <- lapply(step1, numeric)
    step1.character <- lapply(step1, nonnumeric)
    maxelem <- max(sapply(step1, length))
    step1.numeric.t <- lapply(1:maxelem, function(i) sapply(step1.numeric, 
        function(x) x[i]))
    step1.character.t <- lapply(1:maxelem, function(i) sapply(step1.character, 
        function(x) x[i]))
    rank.numeric <- sapply(step1.numeric.t, rank)
    rank.character <- sapply(step1.character.t, 
	function(x) as.numeric(factor(x)))
    rank.numeric[!is.na(rank.character)] <- 0
    rank.character <- t(t(rank.character) + apply(matrix(rank.numeric), 
        2, max, na.rm = TRUE))
    rank.overall <- ifelse(is.na(rank.character), rank.numeric, 
        rank.character)
    order.frame <- as.data.frame(rank.overall)
    if (length(which.nas) > 0) 
        order.frame[which.nas, ] <- Inf
    retval <- do.call("order", order.frame)
    return(retval)
}

writeVECT <- writeVECT6 <- function(SDF, vname, #factor2char = TRUE, 
	v.in.ogr_flags=NULL, ignore.stderr = NULL, driver="ESRI Shapefile") {

        if (get.suppressEchoCmdInFuncOption()) {
            inEchoCmd <- get.echoCmdOption()
             tull <- set.echoCmdOption(FALSE)
        }
        if (is.null(ignore.stderr))
            ignore.stderr <- get("ignore.stderr", envir = .GRASS_CACHE)
        stopifnot(is.logical(ignore.stderr))
    if (requireNamespace("rgdal", quietly = TRUE)) {
        ogrD <- rgdal::ogrDrivers()$name
	if (!(driver %in% ogrD))
            stop(paste("Requested driver", driver, "not available in rgdal"))
        ogrDGRASS <- execGRASS("v.in.ogr", flags="f", intern=TRUE,
            ignore.stderr=ignore.stderr)
        ogrDGRASSs <- strsplit(ogrDGRASS, ": ")
        if (!(driver %in% sapply(ogrDGRASSs, "[", 2)))
            stop(paste("Requested driver", driver, "not available in GRASS"))
        fDrivers <- c("GML", "SQLite")
        dDrivers <- c("ESRI_Shapefile", "MapInfo_File")
        if (!(gsub(" ", "_", driver) %in% c(fDrivers, dDrivers)))
            stop(paste("Requested driver", driver, "not supported"))
        is_dDriver <- TRUE
        if (gsub(" ", "_", driver) %in% fDrivers) is_dDriver <- FALSE
        sss <- strsplit(packageDescription("rgdal")$Version, "-")[[1]]
	type <- NULL
	if (class(SDF) == "SpatialPointsDataFrame") type <- "point"
	if (class(SDF) == "SpatialLinesDataFrame") type <- "line"
	if (class(SDF) == "SpatialPolygonsDataFrame") type <- "boundary"
	if (is.null(type)) stop("Unknown data class")

	pid <- as.integer(round(runif(1, 1, 1000)))
	gtmpfl1 <- dirname(execGRASS("g.tempfile", pid=pid,
            intern=TRUE, ignore.stderr=ignore.stderr))
	rtmpfl1 <- ifelse(.Platform$OS.type == "windows" &&
                (Sys.getenv("OSTYPE") == "cygwin"), 
		system(paste("cygpath -w", gtmpfl1, sep=" "), intern=TRUE), 
		gtmpfl1)

	if (driver == "ESRI Shapefile") 
                shname <- substring(vname, 1, ifelse(nchar(vname) > 8, 8, 
		nchar(vname)))

        if (is_dDriver){
            GDSN <- gtmpfl1
            RDSN <- rtmpfl1
            LAYER <- shname
        } else {
            GDSN <- paste(gtmpfl1, shname, sep=.Platform$file.sep)
            RDSN <- paste(rtmpfl1, shname, sep=.Platform$file.sep)
            LAYER <- shname
        }
        tryCatch(
            {
                rgdal::writeOGR(SDF, dsn=RDSN, layer=LAYER, driver=driver,
                                overwrite_layer=TRUE)
                
                
                execGRASS("v.in.ogr", flags=v.in.ogr_flags,
                          dsn=GDSN, output=vname, layer=LAYER,
                          ignore.stderr=ignore.stderr)
            },
            finally = {
                if (.Platform$OS.type != "windows") {
                    unlink(paste(rtmpfl1, list.files(rtmpfl1, pattern=shname), 
                                 sep=.Platform$file.sep))
                }
            }
        )
        if (get.suppressEchoCmdInFuncOption()) {
            tull <- set.echoCmdOption(inEchoCmd)
        }
    } else {
      stop("rgdal not available")
    }

}

vInfo <- function(vname, layer, ignore.stderr = NULL) {
        if (get.suppressEchoCmdInFuncOption()) {
            inEchoCmd <- get.echoCmdOption()
             tull <- set.echoCmdOption(FALSE)
        }
        if (is.null(ignore.stderr))
            ignore.stderr <- get("ignore.stderr", envir = .GRASS_CACHE)
        stopifnot(is.logical(ignore.stderr))

        if (missing(layer)) layer <- 1L
	vinfo0 <- execGRASS("v.info", flags="t", map=vname,
            layer=layer, intern=TRUE, ignore.stderr=ignore.stderr)

# fix to avoid locale problems 091022

        vinfo1 <- gsub("=", ":", vinfo0)
	con <- textConnection(vinfo1)
	res <- drop(read.dcf(con))
	close(con)
        if (get.suppressEchoCmdInFuncOption()) {
            tull <- set.echoCmdOption(inEchoCmd)
        }
	storage.mode(res) <- "integer"
	res
}

vColumns <- function(vname, layer, ignore.stderr = NULL) {
        if (get.suppressEchoCmdInFuncOption()) {
            inEchoCmd <- get.echoCmdOption()
             tull <- set.echoCmdOption(FALSE)
        }
        if (is.null(ignore.stderr))
            ignore.stderr <- get("ignore.stderr", envir = .GRASS_CACHE)
        stopifnot(is.logical(ignore.stderr))
        if (missing(layer)) layer <- 1L
	vinfo0 <- execGRASS("v.info", flags="c", map=vname,
            layer=layer, intern=TRUE, ignore.stderr=ignore.stderr)       
	con <- textConnection(vinfo0)
        res <- read.table(con, header=FALSE, sep="|")
	close(con)
	names(res) <- c("storageType", "name")
        if (get.suppressEchoCmdInFuncOption()) {
            tull <- set.echoCmdOption(inEchoCmd)
        }
	res
}

vDataCount <- function(vname, layer, ignore.stderr = NULL) {
        if (get.suppressEchoCmdInFuncOption()) {
            inEchoCmd <- get.echoCmdOption()
             tull <- set.echoCmdOption(FALSE)
        }
        if (is.null(ignore.stderr))
            ignore.stderr <- get("ignore.stderr", envir = .GRASS_CACHE)
        stopifnot(is.logical(ignore.stderr))
        column <- "column" %in% parseGRASS("v.db.select")$pnames
        if (missing(layer)) layer <- 1L
        parms <- list(map=vname, layer=layer, columns="cat")
        if (column) tull <- execGRASS("v.db.select", flags="c",
            parameters=parms, intern=TRUE, ignore.stderr=ignore.stderr)
        else tull <- execGRASS("v.db.select", flags="c",
            parameters=parms, intern=TRUE, ignore.stderr=ignore.stderr)
	n <- length(tull)
        if (get.suppressEchoCmdInFuncOption()) {
            tull <- set.echoCmdOption(inEchoCmd)
        }
	n
}


#Date: Thu, 13 Oct 2005 17:34:06 +0200
#From: Markus Neteler <neteler@itc.it>
#
#EXERCISE: HOW LONG ARE COMMON BOUNDARIES OF POLYGONS?
#
#
## Requires: GRASS 6.1-CVS from 13 Oct 2005 or later
##
## data: sudden infant deaths data from North Carolina
## data imported from SHAPE file with v.in.ogr
#
##let's have a look
# d.mon x0
# d.vect sids
#
##we work on a copy:
# g.copy vect=sids,sids_nc
#
##we add a second layer to the map which references the boundaries of
##polygons. In the vector geometry we generate an ID (category) for each
##boundary:
# v.category sids_nc out=sids_nc2 layer=2 type=boundary option=add
#
##Underlying idea:
##we'll fetch the IDs (categories) of the polygons left and right from
##each boundary and store it into the attribute table linked to layer 2.
##In general:
## cat_of_boundary | cat_of_left_polygon | cat_of_right_polygon | length_of_boundary
##
##We want only one category per boundary, that's why the sides check is
##needed (a boundary may consist of several pieces)
##
##So we create a new attribute table and link it to the new layer 2
##of the vector map:
# v.db.addtable sids_nc2 layer=2 col="left integer,right integer,length integer"
#
##Now we query the polygon/boundary relationsships and store it into
##the attribute table linked to layer 2:
# v.to.db map=sids_nc2 option=sides col=left,right layer=2 
#
##Now we have unique categories for the boundaries and can calculate the
##lengths:
# v.to.db map=sids_nc2 option=length col=length layer=2 
#
##Done.
#
##See the new attribute table containing the boundary lengths:
# v.db.select sids_nc2 layer=2
#
## verification (let's check boundary #193):
# d.vect sids_nc2 cat=193 layer=2 col=red type=boundary
# d.zoom
# d.measure
## LEN:     12756.00 meters
#
##what does the attribute table say:
# v.db.select sids_nc2 layer=2 | grep '^193'
##190|65|68|12814
#
##This is reasonably close since on screen digitization in d.measure
##isn't always that precise ...
#

vect2neigh <- function(vname, ID=NULL, ignore.stderr = NULL, remove=TRUE,
    vname2=NULL, units="k") {

    if (get.suppressEchoCmdInFuncOption()) {
        inEchoCmd <- get.echoCmdOption()
        tull <- set.echoCmdOption(FALSE)
    }
    if (is.null(ignore.stderr))
        ignore.stderr <- get("ignore.stderr", envir = .GRASS_CACHE)
    stopifnot(is.logical(ignore.stderr))

    vinfo <- vInfo(vname)
    types <- names(vinfo)[which(vinfo > 0)]
    if (length(grep("areas", types)) == 0) 
		stop("Vector object not of area type")

    n <- vDataCount(vname, ignore.stderr=ignore.stderr)


    if (!is.null(ID)) {
		if (!is.character(ID)) stop("ID not character string")
#		cmd <- paste(paste("v.info", .addexe(), sep=""),
#                    " -c ", vname, sep="")
#		if(.Platform$OS.type == "windows") 
#			tull <- system(cmd, intern=TRUE)
#		else tull <- system(cmd, intern=TRUE, 
#			ignore.stderr=ignore.stderr)
                tull <- execGRASS("v.info", flags="c",
                        map=vname, intern=TRUE, 
			ignore.stderr=ignore.stderr)
		if (length(grep(ID, tull)) == 0)
			stop("ID not found")
#		cmd <- paste(paste("v.db.select", .addexe(), sep=""),
#                    " -c map=", vname, " column=", 
#			ID, sep="")
#		if(.Platform$OS.type == "windows") 
#			ID <- as.character(system(cmd, intern=TRUE))
#		else ID <- as.character(system(cmd, intern=TRUE, 
#			ignore.stderr=ignore.stderr))
                ID <- execGRASS("v.db.select", flags="c", 
                        map=vname, columns=ID, intern=TRUE, 
			ignore.stderr=ignore.stderr) 
		if (length(unique(ID)) != n) 
			stop("fewer than n unique ID values")
    }
    vname2_was_null <- FALSE
    if (is.null(vname2)) {
	
	pid <- as.integer(round(runif(1, 1, 1000)))
	vname2 <- paste(vname, pid, sep="")
#	cmd <- paste(paste("g.copy", .addexe(), sep=""),
#                    " vect=", vname, ",", vname2, sep="")
#	if(.Platform$OS.type == "windows") tull <- system(cmd, intern=TRUE)
#	else tull <- system(cmd, intern=TRUE, ignore.stderr=ignore.stderr)
        tull <- execGRASS("g.copy", vect=paste(vname, 
            vname2, sep=","), intern=TRUE, ignore.stderr=ignore.stderr)
        vname2_was_null <- TRUE
    }
    vname2a <- paste(vname2, "a", sep="")
    if (vname2_was_null) {
#	cmd <- paste(paste("v.category", .addexe(), sep=""),
#                    " ", vname2, " out=", vname2a, 
#		"  layer=2 type=boundary option=add", sep="")
#	if(.Platform$OS.type == "windows") tull <- system(cmd, intern=TRUE)
#	else tull <- system(cmd, intern=TRUE, ignore.stderr=ignore.stderr)
        tull <- execGRASS("v.category", input=vname2,
                output=vname2a, layer=as.integer(2), type="boundary",
                option="add", intern=TRUE, ignore.stderr=ignore.stderr)

#	cmd <- paste(paste("v.db.addtable", .addexe(), sep=""),
#                    " ", vname2a, 
#	" layer=2 col=\"left integer,right integer,length double precision\"", 
#	sep="")
#	if(.Platform$OS.type == "windows") system(cmd)
#	else system(cmd, ignore.stderr=ignore.stderr)
        execGRASS("v.db.addtable", map=vname2a, 
                layer=as.integer(2),
                columns="left integer,right integer,length double precision",
                ignore.stderr=ignore.stderr)

# Using vector map name extended by layer number as table name: landuse175a_2
#Creating table with columns (cat integer, left integer,right integer,length
# double precision)
# The table <landuse175a_2> is now part of vector map <landuse175a> and may
# be deleted or overwritten by GRASS modules
# Select privileges were granted on the table
# Reading features...
# 
# and: no such driver available
# WARNING: Unable to start driver <and>
# ERROR: Unable to open database <C:\Documents> by driver <and>
# Current attribute table links:
# layer <1> table <landuse175a> in database <C:\Documents and Settings\s1155\My Documents\GIS DataBase/Spearfish60/s1155/dbf/> through driver <dbf> with key <cat>
# layer <2> table <landuse175a_2> in database <C:\Documents> through driver <and> with key <cat>
# Vector map <landuse175a@s1155> is connected by:


#	cmd <- paste(paste("v.to.db", .addexe(), sep=""),
#                    " map=", vname2a, 
#		" option=sides col=left,right layer=2", sep="")
#	if(.Platform$OS.type == "windows") system(cmd)
#	else system(cmd, ignore.stderr=ignore.stderr)
        execGRASS("v.to.db", map=vname2a, option="sides",
                columns="left,right", layer=as.integer(2),
                ignore.stderr=ignore.stderr)

#	cmd <- paste(paste("v.to.db", .addexe(), sep=""),
#                    " map=", vname2a, 
#		" option=length col=length layer=2", sep="")
#	if(.Platform$OS.type == "windows") system(cmd)
#	else system(cmd, ignore.stderr=ignore.stderr)
        execGRASS("v.to.db", map=vname2a, option="length",
                columns="length", layer=as.integer(2), units=units,
                ignore.stderr=ignore.stderr)

#	cmd <- paste(paste("v.db.select", .addexe(), sep=""),
#                    " ", vname2a, " layer=2", sep="")
#
#	if(.Platform$OS.type == "windows") res <- system(cmd, intern=TRUE)
#	else res <- system(cmd, intern=TRUE, ignore.stderr=ignore.stderr)
    }
    res <- execGRASS("v.db.select", map=vname2a,
                layer=as.integer(2), intern=TRUE, ignore.stderr=ignore.stderr)

#	cmd <- paste(paste("g.remove", .addexe(), sep=""),
#                    " vect=", vname2, ",", vname2a, sep="")
#	if(.Platform$OS.type == "windows") tull <- system(cmd, intern=TRUE)
#	else tull <- system(cmd, intern=TRUE, ignore.stderr=ignore.stderr)
    if (remove) tull <- execGRASS("g.remove",
            vect=paste(vname2, vname2a, sep=","),
            intern=TRUE, ignore.stderr=ignore.stderr)

    con <- textConnection(res)
    t2 <- read.table(con, sep="|", header=TRUE, row.names=1)
    close(con)
    t3 <- t2[t2$left == -1,]
    t4 <- tapply(t3$length, t3$right, sum)
    external <- numeric(n)
    external[as.integer(names(t4))] <- t4
    t5 <- t2[!t2$left == -1,]
    tmp <- t5$left
    t5$left <- t5$right
    t5$right <- tmp
    t6 <- rbind(t2, t5)
    total <- c(tapply(t6$length, t6$right, sum))
    res <- t6[!t6$left == -1,]
#       avoid integer overflow in by=
#	res <- aggregate(res[3], by=list(left=res$left, right=res$right), sum)
#        dups <- duplicated(res[,1:2])
#        resd <- res[dups,]
#        resda <- aggregate(resd[3], by=list(left=resd$left,
#            right=resd$right), sum)
#        resnd <- res[!dups,]
#        res <- rbind(resda, resnd)
    zz <- paste(res$left, res$right, sep=":")
    uzz <- unique(zz)
    mo <- match(zz, uzz)
    smo <- c(tapply(res[,3], mo, sum))
    names(smo) <- NULL
    suzz <- strsplit(uzz, ":")
    lsuzz <- as.integer(sapply(suzz, "[", 1))
    rsuzz <- as.integer(sapply(suzz, "[", 2))
    reso <- data.frame(left=lsuzz, right=rsuzz, length=smo)
    o <- order(reso$left, reso$right)
    reso <- reso[o,]
    attr(reso, "external") <- external
    attr(reso, "total") <- total
    attr(reso, "region.id") <- ID
    attr(reso, "n") <- n
    class(reso) <- c(class(reso), "GRASSneigh", "spatial.neighbour")
    if (get.suppressEchoCmdInFuncOption()) {
        tull <- set.echoCmdOption(inEchoCmd)
    }

    reso
}

cygwin_clean_temp <- function(verbose=TRUE, ignore.stderr = FALSE) {
	if (Sys.getenv("OSTYPE") != "cygwin") stop("this hack just for cygwin")
	pid <- as.integer(round(runif(1, 1, 1000)))
	cmd <- paste(paste("g.tempfile", .addexe(), sep=""),
                    " pid=", pid, sep="")

	gtmpfl1 <- dirname(ifelse(.Platform$OS.type == "windows", 
		system(cmd, intern=TRUE), system(cmd, 
		intern=TRUE, ignore.stderr=ignore.stderr)))
	rtmpfl1 <- ifelse(.Platform$OS.type == "windows" &&
                (Sys.getenv("OSTYPE") == "cygwin"), 
		system(paste("cygpath -w", gtmpfl1, sep=" "), intern=TRUE), 
		gtmpfl1)
	flst <- list.files(rtmpfl1)
	for (i in flst) {
		f <- paste(rtmpfl1, i, sep="\\")
		x <- file.remove(f)
		if (verbose) cat(f, x, "\n")
	}
	invisible(flst)
}


