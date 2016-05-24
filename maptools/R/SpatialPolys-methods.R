readShapePoly <- function(fn, IDvar=NULL, proj4string=CRS(as.character(NA)), 
	verbose=FALSE, repair=FALSE, force_ring=FALSE, delete_null_obj=FALSE,
	retrieve_ABS_null=FALSE) {
	suppressWarnings(Map <- read.shape(filen=fn, 
		verbose=verbose, repair=repair))
	if (!is.null(IDvar)) {
		IDvar <- as.character(IDvar)
		if (!IDvar %in% names(Map$att.data))
			stop(paste("column not found:", IDvar))
		IDvar <- as.character(Map$att.data[[IDvar]])
	}
	.Map2PolyDF(Map, IDs=IDvar, proj4string=proj4string, 
		force_ring=force_ring, delete_null_obj=delete_null_obj,
		retrieve_ABS_null=retrieve_ABS_null)
}

writePolyShape <- function(x, fn, factor2char = TRUE, max_nchar=254) {
	stopifnot(is(x, "SpatialPolygonsDataFrame"))
	df <- as(x, "data.frame")
	df <- data.frame(SP_ID=I(row.names(df)), df)
	pls <- .SpP2polylist(as(x, "SpatialPolygons"))
	suppressWarnings(write.polylistShape(pls, df, file=fn,
	    factor2char = factor2char, max_nchar=max_nchar))
}

.Map2PolyDF <- function(Map, IDs, proj4string=CRS(as.character(NA)),
	force_ring=FALSE, delete_null_obj=FALSE, retrieve_ABS_null=FALSE) {
# ABS null part shapefiles Graham Williams 080403
# birds NULL part Allen H. Hurlbert 090610
        nullParts <- sapply(Map$Shapes, function(x) x$nParts) == 0
        if (delete_null_obj) {
	    nullParts <- which(nullParts)
	    if (length(nullParts) > 0L) {
              if (!retrieve_ABS_null) {
		for (i in length(nullParts):1)
	            Map$Shapes[[nullParts[i]]] <- NULL
                attr(Map$Shapes,'nshps') <- attr(Map$Shapes,'nshps') - 
                    length(nullParts)
                Map$att.data <- Map$att.data[-nullParts,]
                warning(paste("Null objects with the following", 
                    "indices deleted:", paste(nullParts, collapse=", ")))
              } else {
		res <- Map$att.data[nullParts,]
                return(res)
              }
            }
        } else {
# birds NULL part Allen H. Hurlbert 090610
            if (any(nullParts))
               stop(paste("NULL geometry found:", paste(which(nullParts), collapse=", "), "\n               consider using delete_null_obj=TRUE"))
	}
	if (is.null(IDs))
		IDs <- as.character(sapply(Map$Shapes, function(x) x$shpID))
	SR <- .asSpatialPolygonsShapes(Map$Shapes, IDs, 
		proj4string=proj4string, force_ring=force_ring)
	df <- Map$att.data
	rownames(df) <- IDs
	res <- SpatialPolygonsDataFrame(Sr=SR, data=df)
	res
}

.asSpatialPolygonsShapes <- function(shapes, IDs, 
	proj4string=CRS(as.character(NA)), force_ring=FALSE) {
	if (attr(shapes, "shp.type") != "poly")
		stop("Not polygon shapes")
	if (missing(IDs))
		IDs <- as.character(sapply(shapes, function(x) x$shpID))
	if (length(IDs) != attr(shapes,'nshps')) 
		stop("Number of shapes and IDs differ")
	tab <- table(factor(IDs))
	n <- length(tab)
	IDss <- .mixedsort(names(tab))
# try to preserve sensible ordering
#	IDss <- names(tab)
	reg <- match(IDs, IDss)
	belongs <- lapply(1:n, function(x) which(x == reg))
# assemble the list of Srings
	Srl <- vector(mode="list", length=n)
	for (i in 1:n) {
		nParts <- length(belongs[[i]])
		srl <- NULL
		for (j in 1:nParts) {
			jres <- .shp2srsI(shapes[[belongs[[i]][j]]], 
				.nParts.shpI(shapes[[belongs[[i]][j]]]),
				force_ring=force_ring)
			srl <- c(srl, jres)
		}
		Srl[[i]] <- Polygons(srl, ID=IDss[i])
	}
	res <- as.SpatialPolygons.PolygonsList(Srl, proj4string=proj4string)
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
    if (length(which.blanks) > 0L) 
        x[which.blanks] <- -Inf
    if (length(which.nas) > 0L) 
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
    if (length(which.nas) > 0L) 
        order.frame[which.nas, ] <- Inf
    retval <- do.call("order", order.frame)
    return(retval)
}

.shp2srsI <- function(shp, nParts, force_ring=FALSE) {
	Pstart <- shp$Pstart
	nVerts <- nrow(shp$verts)
	from <- integer(nParts)
	to <- integer(nParts)
	from[1] <- 1
	for (j in 1:nParts) {
		if (j == nParts) to[j] <- nVerts
		else {
			to[j] <- Pstart[j+1]
			from[j+1] <- to[j]+1
		}
	}
	srl <- vector(mode="list", length=nParts)
	for (j in 1:nParts) {
		crds <- shp$verts[from[j]:to[j],,drop=FALSE]
		if (force_ring) {
			if (!isTRUE(identical(crds[1,], crds[nrow(crds),])))
				crds <- rbind(crds, crds[1,])
		}
		srl[[j]] <- Polygon(coords=crds)
	}
	srl
}

.nParts.shpI <- function(shp) attr(shp, "nParts")

.xyList2NAmat <- function(xyList) {
	nParts <- length(xyList)
	res <- xyList[[1]]
	if (nParts > 1) {
		for(i in 2:nParts) 
			res <- rbind(res, c(NA,NA), xyList[[i]])
	}
	res
}

.SpP2polylist <- function(x) {
	pls <- slot(x, "polygons")
	n <- length(pls)
	res <- vector(mode="list", length=n)
	for (i in 1:n) {
		xyL <- lapply(slot(pls[[i]], "Polygons"),
                    function(i) slot(i, "coords"))
		nP <- length(xyL)
		nVs <- sapply(xyL, nrow)
		res[[i]] <- .xyList2NAmat(xyL)
		attr(res[[i]], "nParts") <- as.integer(nP)
		from <- integer(nP)
		to <- integer(nP)
		from[1] <- 1
		to[1] <- nVs[1]
		if (nP > 1) for (j in 2:nP) {
			from[j] <- to[(j-1)] + 2
			to[j] <- from[j] + nVs[j] - 1
		}
		attr(res[[i]], "pstart") <- list(from=as.integer(from), 
			to=as.integer(to))
		attr(res[[i]], "bbox") <- c(bbox(pls[[i]]))
	}
	attr(res, "region.id") <- sapply(pls, function(i) slot(i, "ID"))
	class(res) <- "polylist"
	invisible(res)
}

.polylist2SpP <- function(pl) {
	if (!inherits(pl, "polylist")) stop("not a polylist object")
	n <- length(pl)
	IDs <- attr(pl, "region.id")
	pL <- vector(mode="list", length=n)
	for (i in 1:n) {
		nP <- attr(pl[[i]], "nParts")
		Ps <- vector(mode="list", length=nP)
		from <- attr(pl[[i]], "pstart")$from
		to <- attr(pl[[i]], "pstart")$to
		for (j in 1:nP) {
			Ps[[j]] <- Polygon(pl[[i]][from[j]:to[j],])
		}
		pL[[i]] <- Polygons(Ps, IDs[i])
	}
	res <- SpatialPolygons(pL)
	res
}
