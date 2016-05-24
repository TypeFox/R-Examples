# Copyright (c) 2003-12 by Barry Rowlingson, Roger Bivand, and Edzer Pebesma

getPROJ4VersionInfo <- function() {
    res0 <- .Call("PROJ4VersionInfo", PACKAGE="rgdal")
    res <- paste(res0[[1]], ", [PJ_VERSION: ", res0[[2]], "]", sep="")
    res
}

getPROJ4libPath <- function() {
    res <- Sys.getenv("PROJ_LIB")
    res
}

projNAD <- function() {
    .Call("PROJ4NADsInstalled", PACKAGE="rgdal")
}

"project" <- function(xy, proj, inv=FALSE, use_ob_tran=FALSE) {

    if (!is.numeric(xy)) stop("xy not numeric")
    if (is.matrix(xy)) nc <- dim(xy)[1]
    else if (length(xy) == 2L) nc <- 1
    else stop("xy malformed")
# 111216 RSB
    stopifnot(is.character(proj))
    stopifnot(length(proj) == 1)
    stopifnot(is.logical(inv))
# 120816 RSB
    stopifnot(is.logical(use_ob_tran))
    if (use_ob_tran) {
        gp <- grep("proj=ob_tran", proj)
        if (length(gp) == 0) {
            use_ob_tran <- FALSE
            warning("project: use_ob_tran set FALSE")
# 120820 RSB
        } else inv <- !inv
    }
    if (.Platform$OS.type == "windows" && .Platform$r_arch == "i386") {
     if (!inv) {
        attr(nc, "ob_tran") <- as.integer(use_ob_tran)
        if (attr(nc, "ob_tran") == 0L) {
          res <- .Call("transform", "+proj=longlat", proj, nc,
	    as.double(xy[,1]), as.double(xy[,2]), NULL, PACKAGE="rgdal")
        } else {
          res <- .Call("transform", proj, "+proj=longlat", nc,
	    as.double(xy[,1]), as.double(xy[,2]), NULL, PACKAGE="rgdal")
        }
	if (any(!is.finite(res[[1]])) || any(!is.finite(res[[2]]))) {
	  k <- which(!is.finite(res[[1]]) || !is.finite(res[[2]]))
	  cat("non finite transformation detected:\n")
	  print(cbind(xy, res[[1]], res[[2]])[k,])
	  stop(paste("failure in points", paste(k, collapse=":")))
	}
     } else {
        attr(nc, "ob_tran") <- -as.integer(use_ob_tran)
        if (attr(nc, "ob_tran") == 0L) {
          res <- .Call("transform", proj, "+proj=longlat", nc,
	    as.double(xy[,1]), as.double(xy[,2]), NULL, PACKAGE="rgdal")
        } else {
          res <- .Call("transform", "+proj=longlat", proj, nc,
	    as.double(xy[,1]), as.double(xy[,2]), NULL, PACKAGE="rgdal")
        }
	if (any(!is.finite(res[[1]])) || any(!is.finite(res[[2]]))) {
	  k <- which(!is.finite(res[[1]]) || !is.finite(res[[2]]))
	  cat("non finite transformation detected:\n")
	  print(cbind(xy, res[[1]], res[[2]])[k,])
	  stop(paste("failure in points", paste(k, collapse=":")))
	}
     }
    } else {
     if(!inv) {
# 160404 RSB convert to .Call()
      res <- .Call("project",
                as.integer(nc),
                as.double(xy[,1]),
                as.double(xy[,2]),
                proj,
                as.logical(use_ob_tran),
                PACKAGE="rgdal")
     } else {
      res <- .Call("project_inv",
                as.integer(nc),
                as.double(xy[,1]),
                as.double(xy[,2]),
                proj,
                as.logical(use_ob_tran),
                PACKAGE="rgdal")
     }
    }
    out <- cbind(res[[1]], res[[2]])
    if (!is.null(colnames(xy))) colnames(out) <- colnames(xy)
    out
}


if (!isGeneric("spTransform"))
	setGeneric("spTransform", function(x, CRSobj, ...)
		standardGeneric("spTransform"))

"spTransform.SpatialPoints" <-  function(x, CRSobj, ...) {
	if (is.na(proj4string(x))) 
		stop("No transformation possible from NA reference system")
	if (is.na(CRSargs(CRSobj))) 
		stop("No transformation possible to NA reference system")
	dots = list(...)
        if (!is.null(dots$use_ob_tran)) {
          stopifnot(is.logical(dots$use_ob_tran))
          if (dots$use_ob_tran) {
            gpf <- grep("proj=ob_tran", slot(CRSobj, "projargs"))
            gpi <- grep("proj=ob_tran", proj4string(x))
            if (length(gpf) == 0 && length(gpi) == 0) {
              use_ob_tran <- 0L
              warning("project: use_ob_tran set FALSE")
            } else {
              if (length(gpf) > 0) use_ob_tran <- -1L
              else use_ob_tran <- 1L
            }
          } else {
            use_ob_tran <- 0L
          }
        } else {
          use_ob_tran <- 0L
        }
	crds <- coordinates(x)
	crds.names <- dimnames(crds)[[2]] # crds is matrix
	n <- nrow(crds)
        attr(n, "ob_tran") <- use_ob_tran
        if (ncol(crds) == 2) {
	    res <- .Call("transform", proj4string(x), slot(CRSobj, "projargs"), n,
		as.double(crds[,1]), as.double(crds[,2]), NULL, PACKAGE="rgdal")
	    if (any(!is.finite(res[[1]])) || any(!is.finite(res[[2]]))) {
		k <- which(!is.finite(res[[1]]) || !is.finite(res[[2]]))
		cat("non finite transformation detected:\n")
		print(cbind(crds, res[[1]], res[[2]])[k,])
		stop(paste("failure in points", paste(k, collapse=":")))
	    }
	    crds[,1:2] <- cbind(res[[1]], res[[2]])
        } else {
	    res <- .Call("transform", proj4string(x), slot(CRSobj, "projargs"), n,
		as.double(crds[,1]), as.double(crds[,2]), as.double(crds[,3]),
                PACKAGE="rgdal")
	    if (any(!is.finite(res[[1]])) || any(!is.finite(res[[2]]))
                || any(!is.finite(res[[3]]))) {
		k <- which(!is.finite(res[[1]]) || !is.finite(res[[2]])
                    || !is.finite(res[[3]]))
		cat("non finite transformation detected:\n")
		print(cbind(crds, res[[1]], res[[2]], res[[3]])[k,])
		stop(paste("failure in points", paste(k, collapse=":")))
	    }
	    crds[,1:3] <- cbind(res[[1]], res[[2]], res[[3]])
        }
	# make sure coordinate names are set back:
	dimnames(crds)[[2]] <- crds.names
	x <- SpatialPoints(coords=crds, proj4string=CRSobj)
	x
}
setMethod("spTransform", signature("SpatialPoints", "CRS"), spTransform.SpatialPoints)



"spTransform.SpatialPointsDataFrame" <- function(x, CRSobj, ...) {
	xSP <- as(x, "SpatialPoints")
	resSP <- spTransform(xSP, CRSobj, ...)
	# xDF <- as(x, "data.frame")
	xDF <- x@data # little need to add unique row.names here!
	res <- SpatialPointsDataFrame(coords=coordinates(resSP), data=xDF,
		coords.nrs = numeric(0), proj4string = CRS(proj4string(resSP)))
	res
}
setMethod("spTransform", signature("SpatialPointsDataFrame", "CRS"), 
	spTransform.SpatialPointsDataFrame)




setMethod("spTransform", signature("SpatialPixelsDataFrame", "CRS"), 
	function(x, CRSobj, ...) {
                warning("Grid warping not available, coercing to points")
		spTransform(as(x, "SpatialPointsDataFrame"), CRSobj, ...)})

setMethod("spTransform", signature("SpatialGridDataFrame", "CRS"), 
	function(x, CRSobj, ...) {
                warning("Grid warping not available, coercing to points")
		spTransform(as(x, "SpatialPixelsDataFrame"), CRSobj, ...)})


".spTransform_Line" <- function(x, to_args, from_args, ii, jj,
                use_ob_tran) {
	crds <- slot(x, "coords")
	n <- nrow(crds)
        attr(n, "ob_tran") <- use_ob_tran
	res <- .Call("transform", from_args, to_args, n,
		as.double(crds[,1]), as.double(crds[,2]), NULL,
		PACKAGE="rgdal")
	if (any(!is.finite(res[[1]])) || any(!is.finite(res[[2]]))) {
		k <- which(!is.finite(res[[1]]) || !is.finite(res[[2]]))
		cat("non finite transformation detected:\n")
		print(cbind(crds, res[[1]], res[[2]])[k,])
		stop(paste("failure in Lines", ii, "Line", jj, 
			"points", paste(k, collapse=":")))
	}
	crds <- cbind(res[[1]], res[[2]])
	x <- Line(coords=crds)
	x
}

#setMethod("spTransform", signature("Sline", "CRS"), spTransform.Sline)

".spTransform_Lines" <- function(x, to_args, from_args, ii,
                use_ob_tran) {
	ID <- slot(x, "ID")
	input <- slot(x, "Lines")
	n <- length(input)
	output <- vector(mode="list", length=n)
	for (i in 1:n) output[[i]] <- .spTransform_Line(input[[i]], 
		to_args=to_args, from_args=from_args, ii=ii, jj=i,
                use_ob_tran=use_ob_tran)
	x <- Lines(output, ID)
	x
}

#setMethod("spTransform", signature("Slines", "CRS"), spTransform.Slines)

"spTransform.SpatialLines" <- function(x, CRSobj, ...) {
	from_args <- proj4string(x)
	if (is.na(from_args)) 
		stop("No transformation possible from NA reference system")
	to_args <- slot(CRSobj, "projargs")
	if (is.na(to_args)) 
		stop("No transformation possible to NA reference system")
	dots = list(...)
        if (!is.null(dots$use_ob_tran)) {
          stopifnot(is.logical(dots$use_ob_tran))
          if (dots$use_ob_tran) {
            gpf <- grep("proj=ob_tran", slot(CRSobj, "projargs"))
            gpi <- grep("proj=ob_tran", proj4string(x))
            if (length(gpf) == 0 && length(gpi) == 0) {
              use_ob_tran <- 0L
              warning("project: use_ob_tran set FALSE")
            } else {
              if (length(gpf) > 0) use_ob_tran <- -1L
              else use_ob_tran <- 1L
            }
          } else {
            use_ob_tran <- 0L
          }
        } else {
          use_ob_tran <- 0L
        }
	input <- slot(x, "lines")
	n <- length(input)
	output <- vector(mode="list", length=n)
	for (i in 1:n) output[[i]] <- .spTransform_Lines(input[[i]], 
		to_args=to_args, from_args=from_args, ii=i,
                use_ob_tran=use_ob_tran)
	res <- SpatialLines(output, proj4string=CRS(to_args))
	res
}
setMethod("spTransform", signature("SpatialLines", "CRS"), spTransform.SpatialLines)
"spTransform.SpatialLinesDataFrame" <- function(x, CRSobj, ...) {
	xSP <- as(x, "SpatialLines")
	resSP <- spTransform(xSP, CRSobj, ...)
	xDF <- as(x, "data.frame")
	res <- SpatialLinesDataFrame(sl=resSP, data=xDF, match.ID = FALSE)
	res
}
setMethod("spTransform", signature("SpatialLinesDataFrame", "CRS"), spTransform.SpatialLinesDataFrame)




".spTransform_Polygon" <- function(x, to_args, from_args, ii, jj,
                use_ob_tran) {
	crds <- slot(x, "coords")
	n <- nrow(crds)
        attr(n, "ob_tran") <- use_ob_tran
	res <- .Call("transform", from_args, to_args, n,
		as.double(crds[,1]), as.double(crds[,2]), NULL,
		PACKAGE="rgdal")
	if (any(!is.finite(res[[1]])) || any(!is.finite(res[[2]]))) {
		k <- which(!is.finite(res[[1]]) || !is.finite(res[[2]]))
		cat("non finite transformation detected:\n")
		print(cbind(crds, res[[1]], res[[2]])[k,])
		stop(paste("failure in Polygons", ii, "Polygon", jj, 
			"points", paste(k, collapse=":")))
	}
	crds <- cbind(res[[1]], res[[2]])
	x <- Polygon(coords=crds)
	x
}


".spTransform_Polygons" <- function(x, to_args, from_args, ii,
                use_ob_tran) {
	ID <- slot(x, "ID")
	input <- slot(x, "Polygons")
	n <- length(input)
	output <- vector(mode="list", length=n)
	for (i in 1:n) output[[i]] <- .spTransform_Polygon(input[[i]], 
		to_args=to_args, from_args=from_args, ii=ii, jj=i,
                use_ob_tran=use_ob_tran)
	res <- Polygons(output, ID)
        if (!is.null(comment(x))) comment(res) <- comment(x)
	res
}


"spTransform.SpatialPolygons" <- function(x, CRSobj, ...) {
	from_args <- proj4string(x)
	if (is.na(from_args)) 
		stop("No transformation possible from NA reference system")
	to_args <- slot(CRSobj, "projargs")
	if (is.na(to_args)) 
		stop("No transformation possible to NA reference system")
	dots = list(...)
        if (!is.null(dots$use_ob_tran)) {
          stopifnot(is.logical(dots$use_ob_tran))
          if (dots$use_ob_tran) {
            gpf <- grep("proj=ob_tran", slot(CRSobj, "projargs"))
            gpi <- grep("proj=ob_tran", proj4string(x))
            if (length(gpf) == 0 && length(gpi) == 0) {
              use_ob_tran <- 0L
              warning("project: use_ob_tran set FALSE")
            } else {
              if (length(gpf) > 0) use_ob_tran <- -1L
              else use_ob_tran <- 1L
            }
          } else {
            use_ob_tran <- 0L
          }
        } else {
          use_ob_tran <- 0L
        }
	input <- slot(x, "polygons")
	n <- length(input)
	output <- vector(mode="list", length=n)
	for (i in 1:n) output[[i]] <- .spTransform_Polygons(input[[i]], 
		to_args=to_args, from_args=from_args, ii=i,
                use_ob_tran=use_ob_tran)
	res <- SpatialPolygons(output, pO=slot(x, "plotOrder"), 
		proj4string=CRSobj)
	res
}
setMethod("spTransform", signature("SpatialPolygons", "CRS"), spTransform.SpatialPolygons)

"spTransform.SpatialPolygonsDataFrame" <- function(x, CRSobj, ...) {
	xSP <- as(x, "SpatialPolygons")
	resSP <- spTransform(xSP, CRSobj, ...)
	xDF <- as(x, "data.frame")
	res <- SpatialPolygonsDataFrame(Sr=resSP, data=xDF, match.ID = FALSE)
	res
}
setMethod("spTransform", signature("SpatialPolygonsDataFrame", "CRS"), spTransform.SpatialPolygonsDataFrame)

projInfo <- function(type="proj") {
    opts <- c("proj", "ellps", "datum", "units")
    if (!(type %in% opts)) stop("unknown type")
    t <- as.integer(match(type[1], opts) - 1)
    if (is.na(t)) stop("unknown type")
    res <- .Call("projInfo", t, PACKAGE="rgdal")
    if (type == "proj") res$description <- sapply(strsplit(as.character(
        res$description), "\n"), function(x) x[1])
    res <- data.frame(res)
    res
}
