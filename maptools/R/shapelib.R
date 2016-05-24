# Copyright 2000-2001 (c) Nicholas Lewin-Koh 
# modifications 2001-2008 Roger Bivand
# reads an ESRI shapefile into a map object
# set the variables for the header info


read.shape <- function(filen, dbf.data=TRUE, verbose=TRUE, repair=FALSE) {
  filen <- path.expand(filen)
  .Deprecated("", package="maptools", msg="use readShapeSpatial:\nobjects other than Spatial objects defined in the sp package are deprecated")
  if (length(grep("\\.shp$", tolower(filen))) == 0L)
    filen <- paste(filen, "shp", sep=".")
  shinfo <- getinfo.shape(filen)
  if (dbf.data) {
#    library(foreign)
# filename wrong assumption BDR 100403
#    df <- read.dbf(filen)
    bn <- basename(filen)
    dn <- dirname(filen)
    sbn <- strsplit(bn, "\\.")[[1]]
    lsbn <- length(sbn)
    if (lsbn > 1 && tolower(sbn[lsbn]) == "shp") sbn[lsbn] <- "dbf"
    filen1 <- paste(sbn, collapse=".")
    if (length(grep("\\.dbf$", filen1)) == 0L)
        filen1 <- paste(filen1, "dbf", sep=".")
    if (length(dn) > 0L) {
        filen1 <- paste(dn, filen1, sep=.Platform$file.sep)
    }
    df <- read.dbf(filen1)
    ndf <- as.integer(nrow(df))
  } else ndf <- as.integer(NA)
  if (shinfo[[2]] == 8) {
    if (!dbf.data) stop("to test for multipoint compliance, set dbf.data=TRUE")
    if (ndf != shinfo[[3]]) stop("noncompliant multipoint shapefile")
  }
  shp.lst <- .Call("Rshapeget", as.character(filen), as.logical(repair), 
    PACKAGE="maptools")
  if (verbose) {
    print(shinfo)
  }
  n <- length(shp.lst)
  for (i in 1:n) {
    attr(shp.lst[[i]], "nVerts") <- as.integer(shp.lst[[i]]$nVerts)
    attr(shp.lst[[i]], "nParts") <- as.integer(shp.lst[[i]]$nParts)
    attr(shp.lst[[i]], "shp.type") <- as.integer(shp.lst[[i]]$shp.type)
    attr(shp.lst[[i]], "bbox") <- as.double(shp.lst[[i]]$bbox)
  }
  class(shp.lst) <- "ShapeList"
  if (dbf.data) {
    map <- list(Shapes=shp.lst, att.data=df)
    class(map) <- "Map"
    return(map)
  }
  else {
    return(shp.lst)
  }
}

getinfo.shape <- function(filen) {
  shapehead <-.Call("Rshapeinfo1", as.character(path.expand(filen)), PACKAGE="maptools")
  class(shapehead) <- "shapehead"
  shapehead
}

print.shapehead <- function(x, ...) {
    types <- c("Point", NA, "PolyLine", NA, "Polygon", NA, NA, "MultiPoint", NA, NA, "PointZ", NA, "PolyLineZ", NA, "PolygonZ", NA, NA, "MultiPointZ", NA, NA, "PointM", NA, "PolyLineM", NA, "PolygonM", NA, NA, "MultiPointM", NA, NA, "MultiPatch")
    cat("Shapefile type: ", types[x[[2]]], ", (", x[[2]], "), # of Shapes: ", 
      x[[3]], "\n", sep="")
}


#write.pointShape <- function(object, file, coordinates, factor2char=TRUE, 
write.pointShape <- function(coordinates, df, file, factor2char=TRUE, 
  strictFilename=FALSE, max_nchar=254) {
  .Deprecated("", package="maptools", msg="use writeSpatialShape:\nobjects other than Spatial objects defined in the sp package are deprecated")
  file <- path.expand(file)
  dirnm <- dirname(file)
  bnm0 <- basename(file)
  bnm1 <- strsplit(bnm0, "\\.")[[1]]
  if (bnm1[length(bnm1)] == "shp") 
    bnm <- paste(bnm1[-length(bnm1)], collapse=".")
  else bnm <- bnm0
  file <- paste(dirnm, bnm, sep=.Platform$file.sep)
  if (strictFilename && nchar(basename(file)) > 8) 
    stop("shapefile names must conform to the 8.3 format")
  if (!is.matrix(coordinates)) stop("coordinates must be a matrix")
  if (!is.numeric(coordinates)) stop("coordinates must be numeric")
  ncolcrds <- ncol(coordinates)
  if (ncolcrds < 2) stop("coordinates must have at least 2 columns")
  if (ncolcrds > 3) stop("coordinates must have 2 or 3 columns")
  if (nrow(df) != nrow(coordinates))
    stop("different number of rows in coordinates and data frame")
#  library(foreign)
  write.dbf(df, paste(file, ".dbf", sep=""), factor2char=factor2char, max_nchar=max_nchar)
  storage.mode(coordinates) <- "double"
  res <- .Call("shpwritepoint", as.character(file), coordinates,
    as.integer(ncolcrds), PACKAGE="maptools")
  invisible(res)
}

.isValidPolylist <- function(polylist, verbose=FALSE) {
  if (!inherits(polylist, "polylist")) stop("not a polylist object")
  res <- TRUE
  if (length(polylist) < 1L) {
    if (verbose) cat("zero length polylist\n")
    res <- FALSE
  }
  if (is.null(attr(polylist, "nDims"))) {
    if (verbose) cat("null polylist nDims attribute\n")
    res <- FALSE
  } else {
    if (attr(polylist, "nDims") < 2 || attr(polylist, "nDims") > 3) {
      if (verbose) cat("polylist nDims attribute neither 2 nor 3\n")
      res <- FALSE
    }
    if (!is.integer(attr(polylist, "nDims"))) {
      if (verbose) cat("nDims not all integer\n")
      res <- FALSE
    }
  }
  if (!all(sapply(polylist, function(x) is.double(x)))) {
    if (verbose) cat("coordinates not all double\n")
    res <- FALSE
  }
  if (any(sapply(polylist, function(x) is.null(attr(x, "nParts"))))) {
    if (verbose) cat("null polylist nParts attribute\n")
    res <- FALSE
  } else {
    if (any(sapply(polylist, function(x) attr(x, "nParts") < 1))) {
      if (verbose) cat("polylist nParts attribute less than 1\n")
      res <- FALSE
    }
    if (!all(sapply(polylist, function(x) is.integer(attr(x, "nParts"))))) {
      if (verbose) cat("nParts not all integer\n")
      res <- FALSE
    }
  }
  if (any(sapply(polylist, function(x) is.null(attr(x, "pstart"))))) {
    if (verbose) cat("null polylist pstart attribute\n")
    res <- FALSE
  } else {
    if (any(sapply(polylist, function(x) is.null(attr(x, "pstart")$from)))) {
      if (verbose) cat("null polylist pstart$from attribute\n")
      res <- FALSE
    } else {
      if (!all(sapply(polylist, function(x) is.integer(attr(x, 
        "pstart")$from)))) {
        if (verbose) cat("pstart$from not all integer\n")
        res <- FALSE
      }
    }
    if (any(sapply(polylist, function(x) is.null(attr(x, "pstart")$to)))) {
      if (verbose) cat("null polylist pstart$to attribute\n")
      res <- FALSE
    } else {
      if (!all(sapply(polylist, function(x) is.integer(attr(x, 
        "pstart")$to)))) {
        if (verbose) cat("pstart$to not all integer\n")
        res <- FALSE
      }
    }
  }
  res
}

.makePolylistValid <- function(polylist) {
  if (!inherits(polylist, "polylist")) stop("not a polylist object")
  if (length(polylist) < 1L) stop("zero length polylist")
  n <- length(polylist)
  if (is.null(attr(polylist, "nDims")) || 
    !is.integer(attr(polylist, "nDims")) || 
    (attr(polylist, "nDims") < 2 || attr(polylist, "nDims") > 3)) {
    nD <- unique(sapply(polylist, function(x) dim(x)[2]))
    if (length(nD) > 1L) stop("multiple dimension polylist components")
    nD <- as.integer(nD)
    attr(polylist, "nDims") <- nD
  }
  if (!all(sapply(polylist, function(x) is.double(x)))) {
    for (i in 1:n) { 
      a <- attributes(polylist[[i]])
      polylist[[i]] <- matrix(as.double(polylist[[i]]), ncol=nD)
      attributes(polylist[[i]]) <- a
    }
    warning("coordinates changed to double")
  }
  if (any(sapply(polylist, function(x) is.null(attr(x, "nParts"))))) {
    for (i in 1:n) {
      if (any(is.na(c(polylist[[i]])))) {
	xy <- polylist[[i]]
        NAs <- unclass(attr(na.omit(xy), "na.action"))
	nParts <- length(NAs) + 1L
	from <- integer(nParts)
	to <- integer(nParts)
	from[1] <- 1
	to[nParts] <- nrow(xy)
	if (nParts > 1) {
		for (j in 2:nParts) {
			to[(j-1)] <- NAs[(j-1)]-1
			from[j] <- NAs[(j-1)]+1
		}
	}
        attr(polylist[[i]], "nParts") <- as.integer(nParts)
        a <- list()
	a$from <- as.integer(from)
	a$to <- as.integer(to)
        attr(polylist[[i]], "pstart") <- a
      } else {
        attr(polylist[[i]], "nParts") <- as.integer(1)
        a <- list()
	a$from <- as.integer(1)
	a$to <- as.integer(nrow(polylist[[i]]))
        attr(polylist[[i]], "pstart") <- a
      }
      attr(polylist[[i]], "ringDir") <- as.integer(rep(1,
        attr(polylist[[i]], "nParts")))
      attr(polylist[[i]], "plotOrder") <- 
	as.integer(1:attr(polylist[[i]], "nParts"))
    }
    warning("nParts and pstart added")
  }
  if (any(sapply(polylist, function(x) attr(x, "nParts") < 1)))
    stop("polylist nParts attribute less than 1")
  if (!all(sapply(polylist, function(x) is.integer(attr(x, "nParts"))))) {
    for (i in 1:n) attr(polylist[[i]], "nParts") <- 
		as.integer(attr(polylist[[i]], "nParts"))
    warning("nParts changed to integer")
  }
  if (any(sapply(polylist, function(x) is.null(attr(x, "pstart"))))) {
    for (i in 1:n) {
      if (any(is.na(c(polylist[[i]])))) {
	xy <- polylist[[i]]
        NAs <- unclass(attr(na.omit(xy), "na.action"))
	nParts <- length(NAs) + 1L
	from <- integer(nParts)
	to <- integer(nParts)
	from[1] <- 1
	to[nParts] <- nrow(xy)
	if (nParts > 1) {
		for (j in 2:nParts) {
			to[(j-1)] <- NAs[(j-1)]-1
			from[j] <- NAs[(j-1)]+1
		}
	}
        attr(polylist[[i]], "nParts") <- as.integer(nParts)
        a <- list()
	a$from <- as.integer(from)
	a$to <- as.integer(to)
        attr(polylist[[i]], "pstart") <- a
      } else {
        attr(polylist[[i]], "nParts") <- as.integer(1)
        a <- list()
	a$from <- as.integer(1)
	a$to <- as.integer(nrow(polylist[[i]]))
        attr(polylist[[i]], "pstart") <- a
      }
      attr(polylist[[i]], "ringDir") <- as.integer(rep(1,
        attr(polylist[[i]], "nParts")))
      attr(polylist[[i]], "plotOrder") <- 
	as.integer(1:attr(polylist[[i]], "nParts"))
    }
    warning("nParts and pstart added")
  }
  if (!all(sapply(polylist, function(x) is.integer(attr(x, "pstart")$from)))) {
    for (i in 1:n) attr(polylist[[i]], "pstart")$from <- 
		as.integer(attr(polylist[[i]], "pstart")$from)
    warning("pstart$from changed to integer")
  }
  if (!all(sapply(polylist, function(x) is.integer(attr(x, "pstart")$to)))) {
    for (i in 1:n) attr(polylist[[i]], "pstart")$to <- 
		as.integer(attr(polylist[[i]], "pstart")$to)
    warning("pstart$to changed to integer")
  }
  polylist
}

write.polylistShape <- function(polylist, df, file, factor2char=TRUE, 
  strictFilename=FALSE, force=TRUE, max_nchar=254) {
  .Deprecated("", package="maptools", msg="use writeSpatialShape:\nobjects other than Spatial objects defined in the sp package are deprecated")
  file <- path.expand(file)
  dirnm <- dirname(file)
  bnm0 <- basename(file)
  bnm1 <- strsplit(bnm0, "\\.")[[1]]
  if (bnm1[length(bnm1)] == "shp") 
    bnm <- paste(bnm1[-length(bnm1)], collapse=".")
  else bnm <- bnm0
  file <- paste(dirnm, bnm, sep=.Platform$file.sep)
  if (strictFilename && nchar(basename(file)) > 8) 
    stop("shapefile names must conform to the 8.3 format")
  if (!inherits(polylist, "polylist")) stop("not a polylist object")
  if (length(polylist) < 1L) stop("zero length polylist")
  if (nrow(df) != length(polylist))
    stop("different number of rows in polylist and data frame")
  if (!.isValidPolylist(polylist)) {
    if (!force)
      stop("Invalid polylist - set force=TRUE to coerce to validity")
    else polylist <- .makePolylistValid(polylist)
  }
#  library(foreign)
  write.dbf(df, paste(file, ".dbf", sep=""), factor2char=factor2char, max_nchar=max_nchar)
  res <- .Call("shpwritepolys", as.character(file), polylist, 
    PACKAGE="maptools")
  invisible(res)
}

write.linelistShape <- function(linelist, df, file, factor2char=TRUE, 
  strictFilename=FALSE, max_nchar=254) {
  .Deprecated("", package="maptools", msg="use writeSpatialShape:\nobjects other than Spatial objects defined in the sp package are deprecated")
  file <- path.expand(file)
  dirnm <- dirname(file)
  bnm0 <- basename(file)
  bnm1 <- strsplit(bnm0, "\\.")[[1]]
  if (bnm1[length(bnm1)] == "shp") 
    bnm <- paste(bnm1[-length(bnm1)], collapse=".")
  else bnm <- bnm0
  file <- paste(dirnm, bnm, sep=.Platform$file.sep)
  if (strictFilename && nchar(basename(file)) > 8) 
    stop("shapefile names must conform to the 8.3 format")
  if (length(linelist) < 1L) stop("zero length linelist")
  if (nrow(df) != length(linelist))
    stop("different number of rows in linelist and data frame")
  if (!any(sapply(linelist, function(x) is.integer(attr(x, "nParts"))))) {
    for (i in 1:length(linelist)) { 
      attr(linelist[[i]], "nParts") <- as.integer(attr(linelist[[i]], "nParts"))
    }
    warning("nParts changed to integer")
  }
  if (!any(sapply(linelist, function(x) is.integer(attr(x, "pstart")[[1]])))) {
    for (i in 1:length(linelist)) { 
      attr(linelist[[i]], "pstart")[[1]] <- as.integer(attr(linelist[[i]], 
	"pstart")[[1]])
    }
    warning("pstart changed to integer")
  }
  if (!any(sapply(linelist, function(x) is.integer(attr(x, "pstart")[[2]])))) {
    for (i in 1:length(linelist)) { 
      attr(linelist[[i]], "pstart")[[2]] <- as.integer(attr(linelist[[i]], 
	"pstart")[[2]])
    }
    warning("pstart changed to integer")
  }
#  if (!all(sapply(linelist, function(x) all(!is.na(x)))))
#    stop("NAs in line coordinate data")
  if (!any(sapply(linelist, function(x) is.double(x)))) {
    for (i in 1:length(linelist)) { 
      linelist[[i]] <- matrix(as.double(linelist[[i]]), ncol=2)
    }
    warning("coordinates changed to double")
  }
#  library(foreign)
  write.dbf(df, paste(file, ".dbf", sep=""), factor2char=factor2char, max_nchar=max_nchar)
  res <- .Call("shpwritelines", as.character(file), linelist, 
    PACKAGE="maptools")
  invisible(res)
}

