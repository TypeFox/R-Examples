unionSpatialPolygons <- function(SpP, IDs, threshold=NULL, avoidGEOS=FALSE, avoidUnaryUnion=FALSE) {
    if (!is(SpP, "SpatialPolygons")) stop("not a SpatialPolygons object")
    if (storage.mode(IDs) != "character") IDs <- as.character(IDs)
    if (missing(IDs)) stop("IDs required")
    if (length(slot(SpP, "polygons")) != length(IDs))
        stop("input lengths differ")
    rgeosI <- rgeosStatus()
    if (rgeosI && !avoidGEOS) {
        # require(rgeos)
    	if (!requireNamespace("rgeos", quietly = TRUE))
			stop("package rgeos required for unionSpatialPolygons")
        if (avoidUnaryUnion || rgeos::version_GEOS0() < "3.3.0")
            res <- rgeos::gUnionCascaded(spgeom=SpP, id=IDs)
        else
            res <- rgeos::gUnaryUnion(spgeom=SpP, id=IDs)
    } else {
        stopifnot(isTRUE(gpclibPermitStatus()))
    if (!requireNamespace("gpclib", quietly = TRUE))
		stop("package gpclib required for unionSpatialPolygons")
	# require(gpclib)
	pl <- slot(SpP, "polygons")
	proj4CRS <- CRS(proj4string(SpP))
	SrnParts <- sapply(pl, function(x) length(slot(x, "Polygons")))
	tab <- table(factor(IDs))
	n <- length(tab)
	IDss <- names(tab)
	reg <- match(IDs, IDss)
	belongs <- lapply(1:n, function(x) which(x == reg))
	Srl <- vector(mode="list", length=n)
	for (i in 1:n) {
		ii <- belongs[[i]]
		nParts <- length(ii)
		if (nParts == 1) {
			Srl[[i]] <- Polygons(
				slot(pl[[ii[1]]], "Polygons"), 
				ID=IDss[i])
		} else {
			nPi <- SrnParts[belongs[[i]]]
			m <- sum(nPi)
			pli <- vector(mode="list", length=m)
			jj <- 1
			for (j in 1:nParts) {
				SrSrj <- slot(pl[[ii[j]]], "Polygons")
				for (k in 1:nPi[j]) {
					pli[[jj]] <- slot(SrSrj[[k]], "coords")
					if (jj <= m) jj <- jj + 1
					else stop("jj out of range")
				}
			}
			iin <- length(pli)
			resi <- as(pli[[1]], "gpc.poly")
			for (j in 2:iin) 
				resi <- gpclib::union(resi, as(pli[[j]], "gpc.poly"))
			if (!is.null(threshold)) {
				areas <- sapply(resi@pts, function(x) {
				    gpclib::area.poly(as(cbind(x$x, x$y), "gpc.poly"))})
				resi@pts <- resi@pts[areas > threshold]
			}
			nP <- length(resi@pts)
			Srli <- vector(mode="list", length=nP)
			for (j in 1:nP) {
				crds <- cbind(resi@pts[[j]]$x, resi@pts[[j]]$y)
				crds <- rbind(crds, crds[1,])
				hole <- resi@pts[[j]]$hole
				Srli[[j]] <- Polygon(coords=crds, hole=hole)
			}
			Srl[[i]] <- Polygons(Srli, ID=IDss[i])
		}
	}
	res <- as.SpatialPolygons.PolygonsList(Srl, proj4string=proj4CRS)
    }
    res
}

