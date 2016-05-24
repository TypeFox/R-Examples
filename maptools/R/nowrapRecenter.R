nowrapSpatialPolygons <- function(obj, offset=0, eps=rep(.Machine$double.eps^(1/2.5), 2), avoidGEOS=FALSE) {
    rgeosI <- rgeosStatus()
    if (rgeosI) {
    	if (!requireNamespace("rgeos", quietly = TRUE))
		stop("package rgeos required")
#        require(rgeos)
    } else if (!avoidGEOS) {
    	if (!requireNamespace("polyclip", quietly = TRUE)) {
          warning("package polyclip suggested")
          stopifnot(isTRUE(gpclibPermitStatus()))
    	  if (!requireNamespace("gpclib", quietly = TRUE))
		stop("package gpclib required")
        }
#	require(gpclib)
    }
	if (!is(obj, "SpatialPolygons")) stop("obj not a SpatialPolygons object")
	proj <- is.projected(obj)
	if (is.na(proj)) stop("unknown coordinate reference system")
	if (proj) stop("cannot recenter projected coordinate reference system")
	bblong <- bbox(obj)[1,]
	inout <- bblong[1] < offset && bblong[2] >= offset
	if (inout) {
		pls <- slot(obj, "polygons")
		Srl <- lapply(pls, .nowrapPolygons, offset=offset, eps=eps,
                    rgeosI=rgeosI, avoidGEOS=avoidGEOS)
                Srl <- Srl[!sapply(Srl, is.null)]
		res <- as.SpatialPolygons.PolygonsList(Srl,
			proj4string=CRS(proj4string(obj)))
	} else res <- obj
	res
}

.nowrapPolygons <- function(obj, offset=0, eps=rep(.Machine$double.eps^(1/2.5), 2), rgeosI, avoidGEOS=FALSE) {
	if (!is(obj, "Polygons")) stop("not an Polygons object")
        if (slot(obj, "area")  < mean(eps)) return(NULL)
        if (length(eps) == 1) eps <- rep(eps, 2)
	bbo <- bbox(obj)
	inout <- bbo[1,1] < offset && bbo[1,2] >= offset
	if (inout) {
            if (rgeosI && !avoidGEOS) {
    		if (!requireNamespace("rgeos", quietly = TRUE))
			stop("package rgeos required for .nowrapPolygons")
                comm <- try(rgeos::createPolygonsComment(obj), silent=TRUE)
                isV <- try(rgeos::gIsValid(SpatialPolygons(list(obj))), silent=TRUE)
                 if (class(comm) != "try-error" && class(isV) != "try-error" 
                     && isV) {
                     comment(obj) <- comm
                 } else {
                     stop(paste("invalid Polygons object:", slot(obj, "ID")))
                 }
                 bb <- bbox(obj)
                 bb <- list(x=bb[1,], y=bb[2,])
                 bbmatW <- matrix(c(rep(bb$x[1], 2), rep(offset-eps[1], 2), 
                     bb$x[1], bb$y[1], rep(bb$y[2], 2), rep(bb$y[1], 2)), 
                     ncol=2)
                 bbmatE <- matrix(c(rep(offset+eps[2], 2), rep(bb$x[2], 2), 
                     offset+eps[2], bb$y[1], rep(bb$y[2], 2), 
                     rep(bb$y[1], 2)), ncol=2)
                 SPobj <- SpatialPolygons(list(obj))
                 resW <- rgeos::gIntersection(SPobj, SpatialPolygons(list(Polygons(
                     list(Polygon(bbmatW)), ID="W"))))
                 resE <- rgeos::gIntersection(SPobj, SpatialPolygons(list(Polygons(
                     list(Polygon(bbmatE)), ID="E"))))
                 cparts <- c(slot(slot(resW, "polygons")[[1]], "Polygons"),
                     slot(slot(resE, "polygons")[[1]], "Polygons"))
                 res <- Polygons(cparts, ID=slot(obj, "ID"))
            } else if (isTRUE(gpclibPermitStatus())) {
    		if (!requireNamespace("gpclib", quietly = TRUE))
			stop("package gpclib required for .nowrapPolygons")
		pls <- slot(obj, "Polygons")
		nParts <- length(pls)
		ID <- slot(obj, "ID")
		gpc <- as(slot(pls[[1]], "coords"), "gpc.poly")
		if (nParts > 1) for (i in 2:nParts) gpc <- gpclib::append.poly(gpc, 
			as(slot(pls[[i]], "coords"), "gpc.poly"))
		bb <- gpclib::get.bbox(gpc)
		bbmat1 <- matrix(c(rep(bb$x[1], 2), rep(offset-eps[1], 2), 
			bb$x[1], bb$y[1], rep(bb$y[2], 2), rep(bb$y[1], 2)), 
			ncol=2)
		bbmat2 <- matrix(c(rep(offset+eps[2], 2), rep(bb$x[2], 2), 
			offset+eps[2], bb$y[1], rep(bb$y[2], 2), 
			rep(bb$y[1], 2)), ncol=2)
		gpc_left <- gpclib::intersect(gpc, as(bbmat1, "gpc.poly"))
		gpc_right <- gpclib::intersect(gpc, as(bbmat2, "gpc.poly"))
		gpc_res <- gpclib::append.poly(gpc_left, gpc_right)
		nP <- length(gpc_res@pts)
		if (nP == 0)
			return(obj)
		Srl <- vector(mode="list", length=nP)
		for (j in 1:nP) {
			crds <- cbind(gpc_res@pts[[j]]$x, gpc_res@pts[[j]]$y)
			crds <- rbind(crds, crds[1,])
                        if (any(crds[,1] > offset)) {
                            crds[,1] <- crds[,1] - (2*offset)
                        } else if (any(crds[,1] < -offset)){
                            crds[,1] <- crds[,1] + (2*offset)
                        }
			hole <- gpc_res@pts[[j]]$hole
			rD <- .ringDirxy(crds)
			if (rD == 1 & hole) crds <- crds[nrow(crds):1,]
			if (rD == -1 & !hole)  crds <- crds[nrow(crds):1,]
			Srl[[j]] <- Polygon(coords=crds, hole=hole)
		}
		res <- Polygons(Srl, ID=ID)
            } else {
                if (!requireNamespace("polyclip", quietly = TRUE))
			stop("package polyclip required for .nowrapPolygons")
		pls <- slot(obj, "Polygons")
		nParts <- length(pls)
		ID <- slot(obj, "ID")
                pc <- list(length=nParts)
		crds <- slot(pls[[1]], "coords")
                pc[[1]] <- list(x=crds[,1], y=crds[,2])
		if (nParts > 1) for (i in 2:nParts) {
                        crds <- slot(pls[[i]], "coords")
                        pc[[i]] <- list(x=crds[,1], y=crds[,2])
                }
                bb <- apply(do.call("rbind", lapply(pc, function(i) sapply(i,
                    range))), 2, range)
		bbmat1 <- list(x=unname(c(rep(bb[1, "x"], 2),
                    rep(offset-eps[1], 2), bb[1, "x"])),
                    y=unname(c(bb[1, "y"], rep(bb[2, "y"], 2),
                    rep(bb[1, "y"], 2))))
		bbmat2 <- list(x=unname(c(rep(offset+eps[2], 2),
                    rep(bb[2, "x"], 2), offset+eps[2])),
                    y=unname(c(bb[1, "y"], rep(bb[2, "y"], 2), 
                    rep(bb[1, "y"], 2))))
		pc_left <- polyclip::polyclip(pc, bbmat1, op="intersection")
		pc_right <- polyclip::polyclip(pc, bbmat2, op="intersection")
		pc_res <- c(pc_left, pc_right)
		nP <- length(pc_res)
		if (nP == 0)
			return(obj)
		Srl <- vector(mode="list", length=nP)
		for (j in 1:nP) {
			crds <- cbind(pc_res[[j]]$x, pc_res[[j]]$y)
			crds <- rbind(crds, crds[1,])
                        if (any(crds[,1] > offset)) {
                            crds[,1] <- crds[,1] - (2*offset)
                        } else if (any(crds[,1] < -offset)){
                            crds[,1] <- crds[,1] + (2*offset)
                        }
			Srl[[j]] <- Polygon(coords=crds)
		}
		res <- Polygons(Srl, ID=ID)

            }
	} else res <- obj
	res
}

nowrapRecenter <- function(obj, offset=0, eps=rep(.Machine$double.eps^(1/2.5),
 2), avoidGEOS=FALSE) {
	res <- recenter(nowrapSpatialPolygons(obj, offset=offset, eps=eps,
 avoidGEOS=avoidGEOS))
	res
}


nowrapSpatialLines <- function(obj, offset=0, eps=rep(.Machine$double.eps^(1/2.5), 2)) {
    	if (!requireNamespace("rgeos", quietly = TRUE)) {
	    stop("package rgeos required")
        }

	if (!is(obj, "SpatialLines")) stop("obj not a SpatialLines object")
	proj <- is.projected(obj)
	if (is.na(proj)) stop("unknown coordinate reference system")
	if (proj) stop("cannot recenter projected coordinate reference system")
	proj4CRS <- CRS(proj4string(obj))
	bblong <- bbox(obj)[1,]
	inout <- bblong[1] < offset && bblong[2] >= offset
	if (inout) {
		pls <- slot(obj, "lines")
		Srl <- lapply(pls, .nowrapLines, offset=offset, eps=eps)
		res <- SpatialLines(Srl, proj4CRS)
	} else res <- obj
	res
}


.nowrapLines <- function(obj, offset=0, eps=rep(.Machine$double.eps^(1/2.5), 2)) {
	bbo <- range(sapply(obj@Lines, function(x) range(x@coords[,1])))
	inout <- bbo[1] < offset && bbo[2] >= offset
	if (inout) {
                bb <- bbox(obj)
                bb <- list(x=bb[1,], y=bb[2,])
                bbmatW <- matrix(c(rep(bb$x[1], 2), rep(offset-eps[1], 2), 
                    bb$x[1], bb$y[1], rep(bb$y[2], 2), rep(bb$y[1], 2)), 
                    ncol=2)
                bbmatE <- matrix(c(rep(offset+eps[2], 2), rep(bb$x[2], 2), 
                    offset+eps[2], bb$y[1], rep(bb$y[2], 2), 
                    rep(bb$y[1], 2)), ncol=2)
                SPobj <- SpatialLines(list(obj))
                resW <- rgeos::gIntersection(SPobj, SpatialPolygons(list(
                    Polygons(list(Polygon(bbmatW)), ID="W"))))
                resW <- ifelse (inherits(resW, "SpatialLines"),
                    slot(slot(resW, "lines")[[1]], "Lines"), list(NULL))
                resE <- rgeos::gIntersection(SPobj, SpatialPolygons(list(
                    Polygons(list(Polygon(bbmatE)), ID="E"))))
                resE <- ifelse (inherits(resE, "SpatialLines"),
                    slot(slot(resE, "lines")[[1]], "Lines"), list(NULL))
                cparts <- c(resW, resE)
                cparts <- cparts[!sapply(cparts, is.null)]
                res <- Lines(cparts, ID=slot(obj, "ID"))

	} else res <- obj
	res
}

