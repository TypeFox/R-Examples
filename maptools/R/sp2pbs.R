# PBSmapping utilities

SpatialPolygons2PolySet <- function(SpP) {
	# require(PBSmapping)
    if (!requireNamespace("PBSmapping", quietly = TRUE))
		stop("package PBSmapping required for SpatialPolygons2PolySet")
	pls <- slot(SpP, "polygons")
	n <- length(pls)
	PID <- NULL
	SID <- NULL
	POS <- NULL
	X <- NULL
	Y <- NULL
	for (i in 1:n) {
		srs <- slot(pls[[i]], "Polygons")
		m <- length(srs)
		for (j in 1:m) {
			crds <- slot(srs[[j]], "coords")
			k <- nrow(crds)
			PID <- c(PID, rep(i, k))
			SID <- c(SID, rep(j, k))
#			POS <- c(POS, 1:k)
# Daniel Rodolphe Schlaepfer 140806
			if(slot(srs[[j]], "hole")) POS <- c(POS, k:1) 
			else POS <- c(POS, 1:k)
			X <- c(X, crds[,1])
			Y <- c(Y, crds[,2])
		}
	}
	PID <- as.integer(PID)
	SID <- as.integer(SID)
	POS <- as.integer(POS)
        storage.mode(X) <- "double"
        storage.mode(Y) <- "double"
	pj <- .pbsproj(SpP)
	zn <- NULL
	if (pj == "UTM") {
		zn <- attr(pj, "zone")
		attr(pj, "zone") <- NULL
	}
	res <- PBSmapping::as.PolySet(data.frame(PID=PID, SID=SID, POS=POS, 
		X=X, Y=Y), projection=pj, zone=zn)
	res
}

SpatialLines2PolySet <- function(SL) {
#	require(maps)
	pls <- slot(SL, "lines")
	n <- length(pls)
	PID <- NULL
	SID <- NULL
	POS <- NULL
	X <- NULL
	Y <- NULL
	for (i in 1:n) {
		srs <- slot(pls[[i]], "Lines")
		m <- length(srs)
		for (j in 1:m) {
			crds <- coordinates(srs[[j]])
			k <- nrow(crds)
			PID <- c(PID, rep(i, k))
			SID <- c(SID, rep(j, k))
			POS <- c(POS, 1:k)
			X <- c(X, crds[,1])
			Y <- c(Y, crds[,2])
		}
	}
	PID <- as.integer(PID)
	SID <- as.integer(SID)
	POS <- as.integer(POS)
        storage.mode(X) <- "double"
        storage.mode(Y) <- "double"
	# require(PBSmapping)
    if (!requireNamespace("PBSmapping", quietly = TRUE))
		stop("package PBSmapping required for SpatialPolygons2PolySet")
	pj <- .pbsproj(SL)
	zn <- NULL
	if (pj == "UTM") {
		zn <- attr(pj, "zone")
		attr(pj, "zone") <- NULL
	}
	res <- PBSmapping::as.PolySet(data.frame(PID=PID, SID=SID, POS=POS, 
		X=X, Y=Y), projection=pj, zone=zn)
	res
}

.pbsproj <- function(Sobj) {
	p4str <- proj4string(Sobj)
	if (is.na(p4str)) return("1")
	res <- grep("longlat", p4str, fixed=TRUE)
	if (length(res) > 0L) return("LL")
	res <- regexpr("utm", p4str, fixed=TRUE)
	if (res > 0) {
		val <- "UTM"
		res <- regexpr("+zone=", p4str, fixed=TRUE)
		sres <- substring(p4str, res+attr(res, "match.length"))
		zn0 <- regexpr("[[:digit:]]+", sres)
		attr(val, "zone") <- as.integer(substring(sres, zn0, 
			zn0+attr(zn0, "match.length")))
	} else val <- "1"
	val
}

PolySet2SpatialPolygons <- function(PS, close_polys=TRUE) {
    if (!inherits(PS, "PolySet")) stop("not a PolySet object")
    prj <- attr(PS, "projection")
    if (is.null(prj)) stop("unknown coordinate reference system")
    if (prj == "LL") p4s <- "+proj=longlat +ellps=WGS84"
    else if (prj == "UTM") {
# apparent change in PBS object attributes
        zn <- attr(prj, "zone")
        if (is.null(zn)) zn <- attr(PS, "zone")
        if (is.null(zn)) stop("no valid zone attribute")
	p4s <- paste("+proj=utm +ellps=WGS84 +zone=", zn, sep="")
    } else {
       p4s <- as.character(NA)
       warning("unknown coordinate reference system")
    }
# stop("unknown coordinate reference system") 110310
    hasPID <- "PID" %in% names(PS)
    if (!hasPID) stop("object does not have PID column")
    res0 <- split(PS, PS$PID)
    hasSID <- "SID" %in% names(PS)
    outPolygons <- vector(mode="list", length=length(res0)) 
    if (hasSID) {
        res1 <- lapply(res0, function(x) split(x, x$SID))
        if (close_polys) res1 <- lapply(res1, 
            function(i) lapply(i, function(x) {
                n <- nrow(x)
                if (!isTRUE(identical(x$X[1], x$X[n])) ||
                    !isTRUE(identical(x$Y[1], x$Y[n]))) rbind(x, x[1,])
                else x
            })
        )
# extra level added to fix bug found by A Lobos 080413
        for (i in seq(along=outPolygons)) {
            outPolygons[[i]] <- Polygons(lapply(res1[[i]], function(x) 
                Polygon(cbind(x$X, x$Y))), ID=names(res1)[i])
        }
# PIDs added as IDs 080413
    } else {
        if (close_polys) res0 <- lapply(res0, function(x) {
            n <- nrow(x)
            if (!isTRUE(identical(x$X[1], x$X[n])) ||
                !isTRUE(identical(x$Y[1], x$Y[n]))) rbind(x, x[1,])
            else x
        })
        for (i in seq(along=outPolygons)) {
            outPolygons[[i]] <- Polygons(list(Polygon(cbind(res0[[i]]$X,
                res0[[i]]$Y))), ID=as.character(i))
        }
    }
    outSP <- SpatialPolygons(outPolygons, proj4string=CRS(p4s))
    outSP
}

PolySet2SpatialLines <- function(PS) {
    if (!inherits(PS, "PolySet")) stop("not a PolySet object")
    prj <- attr(PS, "projection")
    prj <- attr(PS, "projection")
    if (is.null(prj)) stop("unknown coordinate reference system")
    if (prj == "LL") p4s <- "+proj=longlat +ellps=WGS84"
    else if (prj == "UTM") {
# apparent change in PBS object attributes
        zn <- attr(prj, "zone")
        if (is.null(zn)) zn <- attr(PS, "zone")
        if (is.null(zn)) stop("no valid zone attribute")
	p4s <- paste("+proj=utm +ellps=WGS84 +zone=", zn, sep="")
    } else {
       p4s <- as.character(NA)
       warning("unknown coordinate reference system")
    }
# stop("unknown coordinate reference system") 110310
    hasPID <- "PID" %in% names(PS)
    if (!hasPID) stop("object does not have PID column")
    res0 <- split(PS, PS$PID)
    hasSID <- "SID" %in% names(PS)
    outLines <- vector(mode="list", length=length(res0)) 
    if (hasSID) {
        res1 <- lapply(res0, function(x) split(x, x$SID))
        for (i in seq(along=outLines)) {
            outLines[[i]] <- Lines(lapply(res1[[i]], function(x) 
                Line(cbind(x$X, x$Y))), ID=as.character(i))
        }
    } else {
        for (i in seq(along=outLines)) {
            outLines[[i]] <- Lines(lapply(res0[[i]], function(x)
                Line(cbind(res0[[i]]$X, res0[[i]]$Y))), ID=as.character(i))
        }
    }
    outSP <- SpatialLines(outLines, proj4string=CRS(p4s))
    outSP
}
