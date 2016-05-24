# Copyright 2003 (c) Barry Rowlingson
# Modified 2006-12 Roger Bivand
###
###
###  Routines for ogr layer data source 
###
###

#
ogrInfo <- function(dsn, layer, encoding=NULL,
  use_iconv=FALSE, swapAxisOrder=FALSE, require_geomType=NULL) {
  if (missing(dsn)) stop("missing dsn")
  if (nchar(dsn) == 0) stop("empty name")
  if (missing(layer)) stop("missing layer")
  if (nchar(layer) == 0) stop("empty name")
  WKB <- c("wkbPoint", "wkbLineString", "wkbPolygon", "wkbMultiPoint",
    "wkbMultiLineString", "wkbMultiPolygon", "wkbGeometryCollection")
  if (!is.null(require_geomType)) {
    stopifnot(is.character(require_geomType) && length(require_geomType)==1)
    m_require_geomType <- match(require_geomType, WKB)
    stopifnot(!is.na(m_require_geomType) || m_require_geomType <= 3)
  }
# a list with various ogr data source information
  
  stopifnot(is.logical(use_iconv))
  stopifnot(length(use_iconv) == 1)
  if (!is.null(encoding)) {
    stopifnot(is.character(encoding))
    stopifnot(length(encoding) == 1)
  }
  if (!use_iconv && !is.null(encoding)) {
    oSE <- getCPLConfigOption("SHAPE_ENCODING")
    tull <- setCPLConfigOption("SHAPE_ENCODING", encoding)
  }
  ogrinfo <- .Call("ogrInfo",as.character(dsn), as.character(layer),
    PACKAGE = "rgdal")
  if (!use_iconv && !is.null(encoding)) {
    tull <- setCPLConfigOption("SHAPE_ENCODING", oSE)
  }

  if (swapAxisOrder) ogrinfo[[5]] <- ogrinfo[[5]][c(2,1,4,3)]

  u_eType <- u_with_z <- null_geometries <- NULL
  deleted_geometries <- NULL
  retain <- NULL
  have_features <- NULL
  all_NULL <- FALSE
  keepGeoms <- NULL

  if (!is.na(ogrinfo[[1]])) {

  fids <- ogrFIDs(dsn=dsn, layer=layer)
  nrows_i <- attr(fids, "i")
  have_features <- nrows_i > 0
  if (have_features && (attr(fids, "i") != attr(fids, "nf"))) {
     retain <- 1:attr(fids, "i")
     afids <- 0:(attr(fids, "nf")-1)
     deleted <- afids[!(afids %in% fids[retain])]
     deleted_geometries <- paste("Deleted feature IDs:", paste(deleted,
        collapse=", "))
     fids <- fids[retain]
  } else {
     deleted_geometries <- NULL
     retain <- NULL
  }
#  attributes(fids) <- NULL
  if (have_features) {
    eTypes <- .Call("R_OGR_types",as.character(dsn), as.character(layer),
      PACKAGE = "rgdal")
    if (is.null(retain)) {
      eType <- eTypes[[4]]
      with_z <- eTypes[[5]]
      isNULL <- as.logical(eTypes[[6]])
   } else {
      eType <- eTypes[[4]][retain]
      with_z <- eTypes[[5]][retain]
      isNULL <- as.logical(eTypes[[6]])[retain]
    }
    null_geometries <- NULL
    if (any(isNULL)) {
      all_NULL <- (sum(isNULL) == length(eType))
      eType <- eType[!isNULL]
      with_z <- with_z[!isNULL]
      null_geometries <- paste("Null geometry IDs:", 
        paste(which(isNULL), collapse=", "))
    }        
    if (!all_NULL) {
      eType[eType == 5L] <- 2L
      eType[eType == 6L] <- 3L

      u_eType <- unique(sort(eType))
      u_with_z <- unique(sort(with_z))
      if (length(u_with_z) != 1L) stop(
        paste("Multiple # dimensions:", 
          paste((u_with_z + 2), collapse=":")))
      if (u_with_z < 0 || u_with_z > 1) stop(
        paste("Invalid # dimensions:", (u_with_z + 2)))

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
          " of ", length(keepGeoms), " features\n",
          "    note that extent applies to all features")
      }
    } else {
      have_features <- FALSE
      warning("ogrInfo: all features NULL")
    }
  }
  } else {
    ogrinfo[[1]] <- attr(ogrinfo[[1]], "dFIDs")
    warning("ogrInfo: feature count overflow")
  }
  names(ogrinfo) <- c("nrows", "nitems", "iteminfo", "driver", "extent",
    "nListFields")
  if (ogrinfo$driver == "ESRI Shapefile") {
      DSN <- dsn
      if (!file.info(DSN)$isdir) DSN <- dirname(normalizePath(dsn))
      DBF_fn <- paste(DSN, .Platform$file.sep, layer, ".dbf", sep = "")
      if (file.exists(DBF_fn)) {
        con <- file(DBF_fn, "rb")
        vr <- readBin(con, "raw", n=32L)
        ldid <- as.integer(vr[30])
        attr(ogrinfo, "LDID") <- ldid
        close(con)
      } else {
        warning("ogrInfo: ", DBF_fn, " not found", sep="")
      }
  }
  names(ogrinfo$iteminfo) <- c("name","type","length","typeName","maxListCount")
  if (use_iconv && !is.null(encoding))
    ogrinfo$iteminfo$name <- iconv(ogrinfo$iteminfo$name, from=encoding)
  ogrinfo$have_features <- have_features
  ogrinfo$eType <- u_eType
  ogrinfo$with_z <- u_with_z
  ogrinfo$null_geometries <- null_geometries
  ogrinfo$deleted_geometries <- deleted_geometries
  ogrinfo$dsn <- dsn
  ogrinfo$layer <- layer
  ogrinfo$p4s <- OGRSpatialRef(dsn, layer)
  if (!is.null(require_geomType))
    attr(ogrinfo, "require_geomType") <- require_geomType
  if (!is.null(keepGeoms)) attr(ogrinfo, "keepGeoms") <- keepGeoms
  class(ogrinfo) <- "ogrinfo"
  ogrinfo
}

print.ogrinfo <- function(x, ...) {
  cat("Source: \"", x$dsn, '\", layer: \"', x$layer, "\"", '\n', sep='')
  cat("Driver:", x$driver)
  if (x$have_features) cat("; number of rows:", x$nrows, "\n")
  WKB <- c("wkbPoint", "wkbLineString", "wkbPolygon", "wkbMultiPoint",
    "wkbMultiLineString", "wkbMultiPolygon", "wkbGeometryCollection")
  if (!is.null(attr(x, "require_geomType"))) {
    cat("  selected geometry type:", attr(x, "require_geomType"), "with",
      sum(attr(x, "keepGeoms")), "rows\n")
  }
  if (!x$have_features) cat(", no features found\n")
  if (x$have_features) cat("Feature type:", paste(WKB[x$eType],
    collapse=", "), "with", x$with_z+2, "dimensions\n")
  if (!is.null(x$extent)) cat("Extent: (", x$extent[1], " ",
    x$extent[2], ") - (", x$extent[3], " ", x$extent[4], ")\n", sep="")
  if (!is.null(x$null_geometries)) cat(x$null_geometries, "\n")
  if (!is.null(x$deleted_geometries)) cat(x$deleted_geometries, "\n")
  if ((nchar(x$p4s) > 1) && !is.na(x$p4s)) cat("CRS:", x$p4s, "\n")
  if (!is.null(attr(x, "LDID"))) cat("LDID:", attr(x, "LDID"), "\n")
  cat("Number of fields:", x$nitems, "\n")
  if (is.null(x$nListFields)) x$nListFields <- 0
  if (x$nListFields > 0) cat("Number of list fields:", x$nListFields, "\n")
  if (x$nitems > 0 && x$nListFields == 0) print(as.data.frame(x$iteminfo)[,1:4])
  if (x$nitems > 0 && x$nListFields > 0 && x$have_features)
    print(as.data.frame(x$iteminfo))
  if (x$nitems > 0 && x$nListFields > 0 && !x$have_features)
    print(as.data.frame(x$iteminfo)[,1:4])
  invisible(x)
}


ogrFIDs <- function(dsn, layer){
  if (missing(dsn)) stop("missing dsn")
  if (nchar(dsn) == 0) stop("empty name")
  if (missing(layer)) stop("missing layer")
  if (nchar(layer) == 0) stop("empty name")
  fids <- .Call("ogrFIDs",as.character(dsn),as.character(layer), PACKAGE = "rgdal")
  if (attr(fids, "i") == 0L) warning("no features found")
  fids
}

ogrDrivers <- function() {
  if (strsplit(getGDALVersionInfo(), " ")[[1]][2] < "2") {
    res <- .Call("ogr_GetDriverNames", PACKAGE="rgdal")
    res <- as.data.frame(res, stringsAsFactors=FALSE)
  } else {
      res <- .Call('RGDAL_GetDriverNames', PACKAGE="rgdal")
      if (!is.null(attr(res, "isVector"))) res$isVector <- attr(res, "isVector")
      res <- as.data.frame(res, stringsAsFactors=FALSE)
      res <- res[res$isVector,]
      names(res)[3] <- "write"
  }
  res <- res[order(res$name),]
  row.names(res) <- NULL
  res
}

"OGRSpatialRef" <- function(dsn, layer) {
    .Call("ogrP4S", as.character(dsn), as.character(layer),
        PACKAGE="rgdal")
}

ogrListLayers <- function(dsn) {
  if (missing(dsn)) stop("missing dsn")
  stopifnot(is.character(dsn))
  stopifnot(length(dsn) == 1)
  if (nchar(dsn) == 0) stop("empty name")
  if (!is.null(attr(dsn, "debug"))) {
    stopifnot(is.logical(attr(dsn, "debug")))
    stopifnot(length(attr(dsn, "debug")) == 1)
  } else {
    attr(dsn, "debug") <- FALSE
  }
  layers <- .Call("ogrListLayers", dsn, PACKAGE = "rgdal")
  n <- length(layers)
  tmp <- layers[n]
  layers <- layers[-n]
  attr(layers, "driver") <- tmp
  attr(layers, "nlayers") <- (n-1)
  layers
}
