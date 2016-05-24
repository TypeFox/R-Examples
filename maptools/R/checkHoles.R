gpclibPermit <- function() {
    if ("gpclib" %in% .packages(all.available = TRUE))
        assign("gpclib", TRUE, envir=.MAPTOOLS_CACHE)
    if (gpclibPermitStatus()) warning("support for gpclib will be withdrawn from maptools at the next major release")
    gpclibPermitStatus()
}

gpclibPermitStatus <- function() get("gpclib", envir=.MAPTOOLS_CACHE)

setRgeosStatus <- function() {
    rgeosI <- "rgeos" %in% .packages(all.available = TRUE)
#    if (rgeosI) {
#        ldNS <- loadedNamespaces()
#        if (!("rgeos" %in% ldNS)) {
#            oo <- try(loadNamespace("rgeos"), silent=TRUE)
#            if (class(oo) == "try-error") rgeosI <- FALSE
#            else unloadNamespace("rgeos")
#        }
#    }
    assign("rgeos", rgeosI, envir=.MAPTOOLS_CACHE)
}

rgeosStatus <- function() get("rgeos", envir=.MAPTOOLS_CACHE)

checkPolygonsHoles <- function(x, properly=TRUE, avoidGEOS=FALSE,
    useSTRtree=FALSE) {
    if (rgeosStatus() && !avoidGEOS) {
        # require(rgeos) # xxx
    	if (!requireNamespace("rgeos", quietly = TRUE))
		stop("package rgeos required")
# version check rgeos
        if (compareVersion(as.character(packageVersion("rgeos")), "0.1-4") < 0)
            useSTRtree <- FALSE
        return(checkPolygonsGEOS(x, properly=properly, useSTRtree=useSTRtree))
    } else {
        stopifnot(isTRUE(gpclibPermitStatus()))
	# require(gpclib)
    	if (!requireNamespace("gpclib", quietly = TRUE))
		stop("package gpclib required")
	if (!is(x, "Polygons")) stop("not an Polygons object")
	pls <- slot(x, "Polygons")
	nParts <- length(pls)
	ID <- slot(x, "ID")
	gpc <- as(slot(pls[[1]], "coords"), "gpc.poly")
	if (nParts > 1) for (i in 2:nParts) gpc <- gpclib::append.poly(gpc, 
		as(slot(pls[[i]], "coords"), "gpc.poly"))
	bb <- gpclib::get.bbox(gpc)
	bbmat <- matrix(c(rep(bb$x[1], 2), rep(bb$x[2], 2), bb$x[1], bb$y[1], 
		rep(bb$y[2], 2), rep(bb$y[1], 2)), ncol=2)
	gpc_bb <- as(bbmat, "gpc.poly")
	gpc_res <- gpclib::intersect(gpc, gpc_bb)
	nP <- length(gpc_res@pts)
	Srl <- vector(mode="list", length=nP)
	for (j in 1:nP) {
		crds <- cbind(gpc_res@pts[[j]]$x, gpc_res@pts[[j]]$y)
		crds <- rbind(crds, crds[1,])
		hole <- gpc_res@pts[[j]]$hole
		rD <- .ringDirxy_gpc(crds)
		if (rD == 1 & hole) crds <- crds[nrow(crds):1,]
		if (rD == -1 & !hole)  crds <- crds[nrow(crds):1,]
		Srl[[j]] <- Polygon(coords=crds, hole=hole)
	}
	res <- Polygons(Srl, ID=ID)
	res
    }
}

.ringDirxy_gpc <- function(xy) {
	a <- xy[,1]
	b <- xy[,2]
	nvx <- length(b)

	if((a[1] == a[nvx]) && (b[1] == b[nvx])) {
		a <- a[-nvx]
		b <- b[-nvx]
		nvx <- nvx - 1
	}
	if (nvx < 3) return(1)

	tX <- 0.0
	dfYMax <- max(b)
	ti <- 1
	for (i in 1:nvx) {
		if (b[i] == dfYMax && a[i] > tX) ti <- i
	}
	if ( (ti > 1) & (ti < nvx) ) { 
		dx0 = a[ti-1] - a[ti]
      		dx1 = a[ti+1] - a[ti]
      		dy0 = b[ti-1] - b[ti]
      		dy1 = b[ti+1] - b[ti]
   	} else if (ti == nvx) {
		dx0 = a[ti-1] - a[ti]
      		dx1 = a[1] - a[ti]
      		dy0 = b[ti-1] - b[ti]
      		dy1 = b[1] - b[ti]
   	} else {
#   /* if the tested vertex is at the origin then continue from 0 (1) */ 
     		dx1 = a[2] - a[1]
      		dx0 = a[nvx] - a[1]
      		dy1 = b[2] - b[1]
      		dy0 = b[nvx] - b[1]
   	}
	v3 = ( (dx0 * dy1) - (dx1 * dy0) )
	if ( v3 > 0 ) return(as.integer(1))
   	else return(as.integer(-1))
}


checkPolygonsGEOS <- function(obj, properly=TRUE, force=TRUE, useSTRtree=FALSE) {
    if (!is(obj, "Polygons")) 
        stop("not a Polygons object")
    if (!requireNamespace("rgeos", quietly = TRUE))
		stop("package rgeos required for checkPolygonsGEOS")
    comm <- try(rgeos::createPolygonsComment(obj), silent=TRUE)
#    isVal <- try(gIsValid(SpatialPolygons(list(obj))), silent=TRUE)
#    if (class(isVal) == "try-error") isVal <- FALSE
    if (class(comm) != "try-error" && !force) {
        comment(obj) <- comm
        return(obj)
    }
    pls <- slot(obj, "Polygons")
    IDs <- slot(obj, "ID")
    n <- length(pls)
    if (n < 1) stop("Polygon list of zero length")
    uniqs <- rep(TRUE, n)
    if (n > 1) {
      if (useSTRtree) tree1 <- rgeos::gUnarySTRtreeQuery(obj)
      SP <- SpatialPolygons(lapply(1:n, function(i) 
        Polygons(list(pls[[i]]), ID=i)))
      for (i in 1:(n-1)) {
        if (useSTRtree) {
            if (!is.null(tree1[[i]])) {
                res <- try(rgeos::gEquals(SP[i,], SP[tree1[[i]],], byid=TRUE),
                    silent=TRUE)
                if (class(res) == "try-error") {
                    warning("Polygons object ", IDs, ", Polygon ",
                        i, ": ", res)
                    next
                }
                if (any(res)) {
                    uniqs[as.integer(rownames(res)[res])] <- FALSE
                }
            }
        } else {
            res <- try(rgeos::gEquals(SP[i,], SP[uniqs,], byid=TRUE), silent=TRUE)
            if (class(res) == "try-error") {
                warning("Polygons object ", IDs, ", Polygon ",
                    i, ": ", res)
                next
            }
            res[i] <- FALSE
            if (any(res)) {
                wres <- which(res)
                uniqs[wres[wres > i]] <- FALSE
            }
        }
      }
    }
    if (any(!uniqs)) warning(paste("Duplicate Polygon objects dropped:",
        paste(wres, collapse=" ")))
    pls <- pls[uniqs]
#    IDs <- IDs[uniqs]
    n <- length(pls)
    if (n < 1) stop("Polygon list of zero length")
    if (n == 1) {
        oobj <- Polygons(pls, ID=IDs)
        comment(oobj) <- rgeos::createPolygonsComment(oobj)
        return(oobj)
    }
    areas <- sapply(pls, slot, "area")
    pls <- pls[order(areas, decreasing=TRUE)]
    oholes <- sapply(pls, function(x) slot(x, "hole"))
    holes <- rep(FALSE, n)
    SP <- SpatialPolygons(lapply(1:n, function(i) 
        Polygons(list(pls[[i]]), ID=i)))
    if (useSTRtree) tree2 <- rgeos::gUnarySTRtreeQuery(SP)
    for (i in 1:(n-1)) {
        if (useSTRtree) {
            if (!is.null(tree2[[i]])) {
                if (properly) res <- rgeos::gContainsProperly(SP[i,], SP[tree2[[i]],],
                    byid=TRUE)
                else res <- rgeos::gContains(SP[i,], SP[tree2[[i]],], byid=TRUE)
            } else {
                res <- FALSE
            }
        } else {
            if (properly) res <- rgeos::gContainsProperly(SP[i,], SP[-(1:i),],
                byid=TRUE)
            else res <- rgeos::gContains(SP[i,], SP[-(1:i),], byid=TRUE)
        }
        wres <- which(res)
        if (length(wres) > 0L) {
            nres <- as.integer(rownames(res))
            holes[nres[wres]] <- ! holes[nres[wres]]
        }
    }
    for (i in 1:n) {
        if (oholes[i] != holes[i])
        pls[[i]] <- Polygon(slot(pls[[i]], "coords"), hole=holes[i])
    }
    oobj <- Polygons(pls, ID=IDs)
    comment(oobj) <- rgeos::createPolygonsComment(oobj)
    oobj    
}
