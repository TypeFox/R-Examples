# Copyright 2006-2012 Roger Bivand

readOGR <- function(dsn, layer, verbose=TRUE, p4s=NULL, 
        stringsAsFactors=default.stringsAsFactors(),
        drop_unsupported_fields=FALSE,
	pointDropZ=FALSE, dropNULLGeometries=TRUE, useC=TRUE,
        disambiguateFIDs=FALSE, addCommentsToPolygons=TRUE, encoding=NULL,
        use_iconv=FALSE, swapAxisOrder=FALSE, require_geomType=NULL,
        integer64="allow.loss") {
	if (missing(dsn)) stop("missing dsn")
	if (nchar(dsn) == 0) stop("empty name")
	if (missing(layer)) stop("missing layer")
	if (nchar(layer) == 0) stop("empty name")
        integer64 <- match.arg(integer64,
          c("allow.loss", "warn.loss", "no.loss"))
        
        int64 <- switch(integer64,
          "allow.loss"=1L,
          "warn.loss"=2L, 
          "no.loss"=3L)
# adding argument for SHAPE_ENCODING environment variable 121124
        stopifnot(is.logical(use_iconv))
        stopifnot(length(use_iconv) == 1)
        if (!is.null(encoding)) {
            stopifnot(is.character(encoding))
            stopifnot(length(encoding) == 1)
        }
        WKB <- c("wkbPoint", "wkbLineString", "wkbPolygon", "wkbMultiPoint",
          "wkbMultiLineString", "wkbMultiPolygon", "wkbGeometryCollection")
        if (!is.null(require_geomType)) {
          stopifnot(is.character(require_geomType) &&
            length(require_geomType)==1)
          m_require_geomType <- match(require_geomType, WKB)
          stopifnot(!is.na(m_require_geomType) || m_require_geomType <= 3)
        }
        
	suppressMessages(ogr_info <- ogrInfo(dsn=dsn, layer=layer,
            encoding=encoding, use_iconv=use_iconv,
            swapAxisOrder=swapAxisOrder, require_geomType=require_geomType))
        if (!ogr_info$have_features) stop("no features found")
        if (is.null(ogr_info$nListFields)) nListFields <- 0
        else nListFields <- ogr_info$nListFields
# 121130 RSB trap no field case (from PostGIS, Mathieu Basille)
        if (ogr_info$nitems > 0) {
          nodata_flag <- FALSE
          if (strsplit(getGDALVersionInfo(), " ")[[1]][2] < "2") {
            keep <- ogr_info$iteminfo$typeName %in% c("Integer", "Real",
              "String", "Date", "Time", "DateTime", "IntegerList",
              "RealList", "StringList")
          } else {
            keep <- ogr_info$iteminfo$typeName %in% c("Integer", "Real",
              "String", "Date", "Time", "DateTime", "IntegerList",
              "RealList", "StringList", "Integer64", "Integer64List")
          }
          if (nListFields > 0)
              ListFields <- as.integer(ogr_info$iteminfo$maxListCount)
          if (drop_unsupported_fields) {
             iflds <- as.integer((1:ogr_info$nitems)-1)
             iflds <- iflds[keep]
             fldnms <- ogr_info$iteminfo$name[keep]
             if (nListFields > 0) ListFields <- ListFields[keep]
             if (any(!keep)) warning(paste("Fields dropped:", 
                 paste(ogr_info$iteminfo$name[!keep], collapse=" ")))
          } else {
             if (any(!keep)) stop(paste("Unsupported field type:", 
                 paste(ogr_info$iteminfo$typeName[!keep], collapse=" ")))
             iflds <- as.integer((1:ogr_info$nitems)-1)
             fldnms <- ogr_info$iteminfo$name
          }
        } else {
          nodata_flag <- TRUE
          iflds <- integer(0)
        }
	fids <- ogrFIDs(dsn=dsn, layer=layer)
        if (attr(fids, "i") != attr(fids, "nf")) {
            retain <- 1:attr(fids, "i")
            afids <- 0:(attr(fids, "nf")-1)
            deleted <- afids[!(afids %in% fids[retain])]
            warning(paste("Deleted feature IDs:", paste(deleted,
                    collapse=", ")))
            fids <- fids[retain]
        } else {
            retain <- NULL
        }
        attributes(fids) <- NULL
# suggestion by Dylan Beaudette 110620
        non_unique_fids <- max(table(fids)) > 1
        if (non_unique_fids) {
            if (disambiguateFIDs) {
                fids <- seq_along(fids) # if present, make new FIDs
            } else {
                stop("FIDs not unique")
            }
        }
	if (verbose) {
		cat("OGR data source with driver:", ogr_info$driver, "\n")
		cat("Source: \"", dsn, '\", layer: \"', layer, "\"", '\n',
			sep='')
		cat("with", length(fids), "features")
                if (!is.null(attr(ogr_info, "require_geomType")))
                    cat(";\nSelected", attr(ogr_info, "require_geomType"),
                      "feature type, with", sum(attr(ogr_info, "keepGeoms")),
                      "rows")
                cat("\n")
		cat("It has", length(iflds), "fields")
                if (nListFields > 0)
                  cat(", of which", nListFields, "list fields")
                cat("\n")
	}
# suggestion by Paul Hiemstra 070817
	prj <- .Call("ogrP4S", as.character(dsn), as.character(layer), 
		PACKAGE="rgdal")
	if (!is.null(p4s)) {
            if (!is.na(prj)) {
                warning("p4s= argument given as: ", p4s, "\n and read as: ", prj, 
                "\n read string overridden by given p4s= argument value")
            }
        } else {
            p4s <- prj
        }

	if (!is.na(p4s) && nchar(p4s) == 0) p4s <- as.character(NA)

# adding argument for SHAPE_ENCODING environment variable 121124
        if (!use_iconv && !is.null(encoding) && 
            ogr_info$driver == "ESRI Shapefile") {
            oSE <- getCPLConfigOption("SHAPE_ENCODING")
            tull <- setCPLConfigOption("SHAPE_ENCODING", encoding)
        }
	if (nodata_flag) {
            dlist <- list(FID=as.integer(fids))
        } else {
            attr(iflds, "nListFields") <- as.integer(nListFields)
            nflds <- length(iflds)
            if (nListFields > 0) {
                attr(iflds, "ListFields") <- ListFields
                nflds <- nflds + sum(ListFields) - nListFields
                fldnms1 <- NULL
                for (i in seq(along=ListFields)) {
                    if (ListFields[i] == 0) fldnms1 <- c(fldnms1, fldnms[i])
                    else fldnms1 <- c(fldnms1,
                        paste(fldnms[i], 1:ListFields[i], sep=""))
                }
                stopifnot(length(fldnms1) == nflds)
                fldnms <- fldnms1
            }
            attr(iflds, "nflds") <- as.integer(nflds)
            attr(iflds, "int64") <- as.integer(int64)
            dlist <- .Call("ogrDataFrame", as.character(dsn),
                as.character(layer), as.integer(fids), iflds, PACKAGE="rgdal")
	    names(dlist) <- make.names(fldnms ,unique=TRUE)

            if (use_iconv && !is.null(encoding)) {
                for (i in seq(along=dlist)) {
                    if (is.character(dlist[[i]])) {
                       dlist[[i]] <- iconv(dlist[[i]], from=encoding)
                    }
                }
            }
        }
        if (!use_iconv && !is.null(encoding) && 
            ogr_info$driver == "ESRI Shapefile") {
            tull <- setCPLConfigOption("SHAPE_ENCODING", oSE)
        }


	geometry <- .Call("R_OGR_CAPI_features", as.character(dsn), 
		as.character(layer), comments=addCommentsToPolygons,
                PACKAGE="rgdal")
	if (is.null(retain)) {
	    eType <- geometry[[4]]
	    with_z <- geometry[[6]]
            isNULL <- as.logical(geometry[[7]])
	    gFeatures <- geometry[[5]]
            gComments <- geometry[[8]]
        } else {
	    eType <- geometry[[4]][retain]
	    with_z <- geometry[[6]][retain]
            isNULL <- as.logical(geometry[[7]])[retain]
	    gFeatures <- geometry[[5]][retain]
            gComments <- geometry[[8]][retain]
        }
        rm(geometry);
        gc(verbose = FALSE)
        if (any(isNULL)) {
            eType <- eType[!isNULL]
            with_z <- with_z[!isNULL]
        }   
     
	u_with_z <- unique(sort(with_z))
	if (length(u_with_z) != 1L) stop(
		paste("Multiple # dimensions:", 
			paste((u_with_z + 2), collapse=":")))
	if (u_with_z < 0 || u_with_z > 1) stop(
		paste("Invalid # dimensions:", (u_with_z + 2)))

        eType[eType == 5L] <- 2L
        eType[eType == 6L] <- 3L

	u_eType <- unique(sort(eType))

    t_eType <- table(eType)
    if (is.null(require_geomType)) {
      keepGeoms <- NULL
      if (length(u_eType) > 1L) stop(
        paste("Multiple incompatible geometries:", 
          paste(paste(WKB[as.integer(names(t_eType))], t_eType, sep=": "),
          collapse="; ")))
#  if (length(u_eType) == 2L) {
#    if (u_eType[1] == 2 && u_eType[2] == 5) u_eType = 2
#    else if (u_eType[1] == 3 && u_eType[2] == 6) u_eType = 3
#    else stop(paste("Multiple incompatible geometries:", 
#      paste(paste(WKB[as.integer(names(t_eType))], t_eType, sep=": "),
#        collapse="; ")))
#   }
    } else {
      if (!require_geomType %in% WKB[as.integer(names(t_eType))])
        stop(require_geomType, "not in", WKB[as.integer(names(t_eType))])
      u_eType <- match(require_geomType, WKB)
      keepGeoms <- WKB[eType] == require_geomType
      message("NOTE: keeping only ", sum(keepGeoms), " ", require_geomType,
        " of ", length(keepGeoms), " features\n")
    }

	data <- data.frame(dlist, row.names=fids,
            stringsAsFactors=stringsAsFactors)
        rm(dlist)
        gc(verbose = FALSE)
	if (length(gFeatures) != length(fids)) stop("Feature mismatch")

        if (any(isNULL)) {
            if (dropNULLGeometries) {
                warning(paste("Dropping null geometries:", paste(which(isNULL),
                    collapse=", ")))
                gFeatures <- gFeatures[!isNULL]
	        data <- data[!isNULL, , drop=FALSE]
                fids <- fids[!isNULL]
                gComments <- gComments[!isNULL]
            } else {
                warning(paste("Null geometries found:", paste(which(isNULL),
                    collapse=", ")))
                warning("dropNULLGeometries FALSE, returning only data for null-geometry features")
                return(data[isNULL, , drop=FALSE])
            }
        }

        if (!is.null(require_geomType)) {
                gFeatures <- gFeatures[keepGeoms]
	        data <- data[keepGeoms, , drop=FALSE]
                fids <- fids[keepGeoms]
                gComments <- gComments[keepGeoms]
        }

	if (u_eType == 1) { # points
		if (u_with_z == 0 || pointDropZ) {
                    if (swapAxisOrder) {
			coords <- do.call("rbind", lapply(gFeatures, 
				function(x) c(x[[1]][[2]], x[[1]][[1]])))
                    } else {
			coords <- do.call("rbind", lapply(gFeatures, 
				function(x) c(x[[1]][[1]], x[[1]][[2]])))
                    }
		} else {
                    if (swapAxisOrder) {
			coords <- do.call("rbind", lapply(gFeatures, 
				function(x) c(x[[1]][[2]], x[[1]][[1]],
				x[[1]][[3]])))
                    } else {
			coords <- do.call("rbind", lapply(gFeatures, 
				function(x) c(x[[1]][[1]], x[[1]][[2]],
				x[[1]][[3]])))
                    }  
		}
#		data <- data.frame(dlist)
		row.names(data) <- NULL
		res <- SpatialPointsDataFrame(coords=coords, data=data,
			proj4string=CRS(p4s))
	} else if (u_eType == 2) { # lines
		if (u_with_z != 0) warning("Z-dimension discarded")
		n <- length(gFeatures)
		lnList <- vector(mode="list", length=n)
		for (i in 1:n) {
			iG <- gFeatures[[i]]
			m <- length(iG)
			lnlist <- vector(mode="list", length=m)
			for (j in 1:m) {
				jG <- iG[[j]]
                                if (swapAxisOrder) {
				  lnlist[[j]] <- Line(cbind(jG[[2]], jG[[1]]))
                                } else {
				  lnlist[[j]] <- Line(cbind(jG[[1]], jG[[2]]))
                                }
			}
			lnList[[i]] <- Lines(lnlist, ID=as.character(fids[i]))
		}
		SL <- SpatialLines(lnList, proj4string=CRS(p4s))
#		data <- data.frame(dlist, row.names=fids)
		res <- SpatialLinesDataFrame(SL, data)
	} else if (u_eType == 3) { # polygons
            if (u_with_z != 0) warning("Z-dimension discarded")
            if (useC) {
#                plList <- .Call("make_polygonslist", gFeatures,
#                    as.list(as.character(fids)), PACKAGE="rgdal")
		n <- length(gFeatures)
		plList <- vector(mode="list", length=n)
		for (i in 1:n) {
			iG <- gFeatures[[i]]
                        if (swapAxisOrder) {
                          iG <- lapply(iG, function(x) {
                            tmp <- x[[1]]; x[[1]] <- x[[2]]; x[[2]] <- tmp
                          })
                        }
                        if (addCommentsToPolygons) {
                            thisPL <- Polygons(.Call("make_Polygonlist",
                                iG, gComments[[i]], PACKAGE="rgdal"),
                                ID=as.character(fids[i]))
                            comment(thisPL) <- paste(gComments[[i]],
                                collapse=" ")
                        } else {
                            thisPL <- Polygons(.Call("make_Polygonlist",
                                iG, NULL, PACKAGE="rgdal"),
                                ID=as.character(fids[i]))
                        }
			plList[[i]] <- thisPL
                    }
            } else {
		n <- length(gFeatures)
		plList <- vector(mode="list", length=n)
		for (i in 1:n) {
			iG <- gFeatures[[i]]
                        if (swapAxisOrder) {
                          iG <- lapply(iG, function(x) {
                            tmp <- x[[1]]; x[[1]] <- x[[2]]; x[[2]] <- tmp
                          })
                        }
			m <- length(iG)
			pllist <- vector(mode="list", length=m)
			for (j in 1:m) {
				jG <- iG[[j]]
				cmat <- cbind(jG[[1]], jG[[2]])
				if (!identical(cmat[1,], cmat[nrow(cmat),])) {
				  cmat <- rbind(cmat, cmat[1,])
                                  warning(paste("Ring closed in Polygons",
				    i, "Polygon", j))
				}
				t0 <- try(pllist[[j]] <- Polygon(cmat),
                                    silent=TRUE)
                                if (class(t0) == "try-error") {
                                     print(cmat)
                                     print(t0)
                                     stop("i: ", i, ", j: ", j,
                                       ", Polygon error exit")
                                }
			}
                        thisPL <- Polygons(pllist, ID=as.character(fids[i]))
                        if (addCommentsToPolygons) {
                            comment(thisPL) <- paste(gComments[[i]],
                                collapse=" ")
                            if (!isTRUE(all.equal(as.logical(gComments[[i]]),
                              sapply(slot(thisPL, "Polygons"), slot, "hole"))))
                              warning("comment/hole mismatch, geometry:", i)
                        }
			plList[[i]] <- thisPL
                    }
		}
                rm(gFeatures)
                gc(verbose = FALSE)
		SP <- SpatialPolygons(plList, proj4string=CRS(p4s))
                rm(plList)
                gc(verbose = FALSE)
#		data <- data.frame(dlist, row.names=fids)
		res <- SpatialPolygonsDataFrame(SP, data, match.ID=FALSE)
	} else stop(paste("Incompatible geometry:", u_eType))

	res
}

showWKT <- function(p4s, file=NULL, morphToESRI=TRUE) {

	if (!is.character(p4s)) stop("invalid p4s object")
        stopifnot(length(p4s) == 1)
	if (!is.logical(morphToESRI)) stop("invalid morphToESRI object")
	res <- .Call("p4s_to_wkt", as.character(p4s), as.integer(morphToESRI), 
		PACKAGE="rgdal")
	if (!is.null(file)) cat(res, "\n", sep="", file=file)
	res
}

showP4 <- function(wkt, morphFromESRI=TRUE) {

	if (!is.character(wkt)) stop("invalid wkt object")
        stopifnot(length(wkt) == 1)
	if (!is.logical(morphFromESRI)) stop("invalid morphFromESRI object")
	res <- .Call("wkt_to_p4s", as.character(wkt),
                as.integer(morphFromESRI), PACKAGE="rgdal")
	res
}


showEPSG <- function(p4s) {

	if (!is.character(p4s)) stop("invalid p4s object")
	res <- .Call("ogrAutoIdentifyEPSG", as.character(p4s), PACKAGE="rgdal")
	res
}

 
