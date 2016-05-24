writeOGR <- function(obj, dsn, layer, driver, dataset_options=NULL, layer_options=NULL, verbose=FALSE, check_exists=NULL, overwrite_layer=FALSE, delete_dsn=FALSE, morphToESRI=NULL, encoding=NULL) {
    stopifnot(is.logical(verbose))
    drvs <- ogrDrivers()
    mch <- match(driver, drvs$name)
    if (is.na(mch) || length(mch) > 1L)
        stop(paste("No such driver:", driver))
    if (!drvs$write[mch]) stop("Chosen driver cannot create files")
    if (driver == "KML" || driver == "GPX") {
      if (is.na(is.projected(obj))) warning(paste("Unknown coordinate", 
        "reference system: for KML driver should be geographical"))
      else if (is.projected(obj)) {
# added brackets 120330
         warning(paste("Projected coordinate",
        "reference system: for KML driver should be geographical"))
        proj4string(obj) <- CRS(as.character(NA))
      }
# fix for over-eager internal checking in the KML driver 090114
    }
    stopifnot(inherits(obj, "Spatial"))
    if (gridded(obj)) {
        obj <- as(obj, "SpatialPointsDataFrame")
        if (verbose) warning("coercing gridded object to points")
    }
    if (!"data" %in% names(getSlots(class(obj))))
        stop("obj must be a SpatialPointsDataFrame, SpatialLinesDataFrame or\n    SpatialPolygonsDataFrame") 
    dfcls <- sapply(slot(obj, "data"), function(x) class(x)[1])
    dfcls <- gsub("POSIXct", "POSIXt", dfcls)
    dfcls <- gsub("POSIXlt", "POSIXt", dfcls)
    dfcls <- gsub("Date", "POSIXt", dfcls)
# fix for logical and better reporting Barry Rowlingson 091106
    known <- c("numeric", "character", "factor", "POSIXt", "integer", "logical")
    if (!all(dfcls %in% known)) 
        stop("Can't convert columns of class: ",
            paste(unique(dfcls[!dfcls %in% known]), collapse=","),
            "; column names: ", paste(names(obj@data)[!dfcls %in% known],
            collapse=","))
    dftof <- sapply(slot(obj, "data"), typeof)
    known <- c("double", "character", "integer", "logical", "list")
    if (!all(dftof %in% known))
        stop("Can't convert columns of type: ", 
            paste(unique(dftof[!dftof %in% known]),collapse=","), 
            "; column names: ", paste(names(obj@data)[!dftof %in% known], 
            collapse=","))
    if (any("list" %in% dftof)) {
         whch <- dfcls[which(dftof == "list")] != "POSIXt"
         if (any(whch)) 
             stop("Can't convert columns of type list:",
                 paste(names(obj@data)[whch], collapse=","))
    }

    if (is.null(check_exists)) {
        if (as.integer(getGDALVersionInfo("VERSION_NUM")) >= 1800L) {
            check_exists <- TRUE
        } else {
            check_exists <- FALSE
        }
    }

    odsn <- NULL
    if (check_exists) {
        already_exists <- FALSE
        ogrI <- .Call("ogrCheckExists", as.character(dsn),
            as.character(layer), PACKAGE = "rgdal")
        if (ogrI) {
            if (driver == attr(ogrI, "driver")) already_exists <- TRUE
        }
        if (already_exists) {
            if (overwrite_layer) {
# NASTY KLUDGE: GDAL 2 appears to handle shapefile deletion awkwardly
# and prefers the "single" *.shp definition of dsn
                if (strsplit(getGDALVersionInfo(), " ")[[1]][2] >= "2" && 
                    driver == "ESRI Shapefile" && file.info(dsn)$isdir) {
                    odsn <- dsn
                    dsn <- paste(dsn, "/", layer, ".shp", sep="")
                    if (verbose) warning(dsn, " substituted for ", odsn,
                        " for layer deletion")
                    odelete_dsn <- delete_dsn
                    delete_dsn <- TRUE
                    owarnZZ <- options("warn")$warn
                    options(warn=-1)
                }
                layer_del <- try(ogrDeleteLayer(dsn, layer, driver),
                    silent=TRUE)
                if (class(layer_del) == "try-error" && delete_dsn) {
                    ogrDeleteDataSource(dsn, driver)
                    if (verbose) warning("existing data source removed")
                }
                if (exists("odelete_dsn")) {
                    delete_dsn <- odelete_dsn
                    options(warn=owarnZZ)
                }
                if (verbose) warning("existing layer removed")
            } else {
                stop("layer exists, use a new layer name")
            }
        }
    }
    if (is.null(morphToESRI))
        morphToESRI <- ifelse(driver == "ESRI Shapefile", TRUE, FALSE)
    stopifnot(is.logical(morphToESRI))
    stopifnot(length(morphToESRI) == 1)

    nf <- length(dfcls)
    ldata <- vector(mode="list", length=nf)
    ogr_ftype <- integer(nf)
    for (i in 1:nf) {
        if (dfcls[i] == "numeric" && dftof[i] == "double") {
            ldata[[i]] <- slot(obj, "data")[,i]
            ogr_ftype[i] <- as.integer(2) #"OFTReal"
        } else if (dfcls[i] == "character" && dftof[i] == "character") {
            ldata[[i]] <- slot(obj, "data")[,i]
            ogr_ftype[i] <- as.integer(4) #"OFTString"
        } else if (dfcls[i] == "factor" && dftof[i] == "integer") {
            ldata[[i]] <- as.character(slot(obj, "data")[,i])
            ogr_ftype[i] <- as.integer(4) #"OFTString"
        } else if (dfcls[i] == "POSIXt" && dftof[i] == "integer") {
            ldata[[i]] <- as.character(format(slot(obj, "data")[,i]))
            ogr_ftype[i] <- as.integer(4) #"OFTString"
        } else if (dfcls[i] == "POSIXt" && dftof[i] == "double") {
            ldata[[i]] <- as.character(format(slot(obj, "data")[,i]))
            ogr_ftype[i] <- as.integer(4) #"OFTString"
        } else if (dfcls[i] == "POSIXt" && dftof[i] == "list") {
            ldata[[i]] <- as.character(format(slot(obj, "data")[,i]))
            ogr_ftype[i] <- as.integer(4) #"OFTString"
        } else if (dfcls[i] == "integer" && dftof[i] == "integer") {
            ldata[[i]] <- slot(obj, "data")[,i]
            ogr_ftype[i] <- as.integer(0) #"OFTInteger"
        } else if (dfcls[i] == "logical" && dftof[i] == "logical") {
# fix for logical and better reporting Barry Rowlingson 091106
            ldata[[i]] <- as.integer(slot(obj, "data")[,i])
            ogr_ftype[i] <- as.integer(0) #"OFTInteger"
        } else stop(paste(dfcls[i], dftof[i], "unknown data type"))
        if (!is.null(encoding)) {
            if (ogr_ftype[i] == 4L) {
                ldata[[i]] <- iconv(ldata[[i]], from=encoding, to="UTF-8")
            }
        }
    }
    fld_names <- names(dfcls)
    if (!is.null(encoding)) {
        fld_names <- iconv(fld_names, from=encoding, to="UTF-8")
    }
    if (driver == "ESRI Shapefile") {
        if (any(nchar(fld_names) > 10)) {
            fld_names <- abbreviate(fld_names, minlength=7)
            warning("Field names abbreviated for ESRI Shapefile driver")
            if (any(nchar(fld_names) > 10)) 
                fld_names <- abbreviate(fld_names, minlength=5)
        }
# fix for dots in DBF field names 121124
        if (length(wh. <- grep("\\.", fld_names) > 0)) {
            fld_names[wh.] <- gsub("\\.", "_", fld_names[wh.])
        }
    }
    if (length(fld_names) != length(unique(fld_names)))
       stop("Non-unique field names")
    nobj <- nrow(slot(obj, "data"))
# add FIDs 130502
    owarn <- options("warn")
    options("warn"=-1L)
    FIDs <- as.integer(row.names(obj))
    if (any(is.na(FIDs))) FIDs <- as.integer(0:(nobj-1))
    options("warn"=owarn$warn)
    attr(nf, "verbose") <- as.logical(verbose)
    
    pre <- list(obj, as.character(dsn), as.character(layer), 
        as.character(driver), as.integer(nobj), nf,
        as.character(fld_names), as.integer(ogr_ftype), ldata, 
        as.character(dataset_options), as.character(layer_options),
        as.logical(morphToESRI), as.integer(FIDs))
    res <- .Call("OGR_write", pre, PACKAGE="rgdal")
    if (verbose) {
        res <- list(object_type=res, output_dsn=dsn, output_layer=layer,
            output_diver=driver, output_n=nobj, output_nfields=nf,
            output_fields=fld_names, output_fclasses=ogr_ftype, 
            dataset_options=dataset_options, layer_options=layer_options,
            morphToESRI=morphToESRI, FIDs=FIDs)
        return(res)
# add FIDs 130502
    } 
}

ogrDeleteLayer <- function(dsn, layer, driver) {
    if (missing(dsn)) stop("missing dsn")
    if (nchar(dsn) == 0) stop("empty name")
    if (missing(layer)) stop("missing layer")
    if (nchar(layer) == 0) stop("empty name")
    if (missing(driver)) stop("missing layer")
    if (nchar(driver) == 0) stop("empty name")
    res <- .Call("ogrDeleteLayer", as.character(dsn), as.character(layer),
        as.character(driver), PACKAGE = "rgdal")
    invisible(res)
}

ogrDeleteDataSource <- function(dsn, driver) {
    if (missing(dsn)) stop("missing dsn")
    if (nchar(dsn) == 0) stop("empty name")
    if (missing(driver)) stop("missing layer")
    if (nchar(driver) == 0) stop("empty name")
    res <- .Call("ogrDeleteDataSource", as.character(dsn), as.character(driver),
        PACKAGE = "rgdal")
    invisible(res)

}

