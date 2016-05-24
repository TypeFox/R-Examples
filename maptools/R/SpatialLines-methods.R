readShapeLines <- function(fn, proj4string=CRS(as.character(NA)), 
	verbose=FALSE, repair=FALSE, delete_null_obj=FALSE) {
	suppressWarnings(Map <- read.shape(filen=fn, verbose=verbose,
	    repair=repair))
	suppressWarnings(.shp2LinesDF(Map, proj4string=proj4string,
            delete_null_obj=delete_null_obj))
}

writeLinesShape <- function(x, fn, factor2char = TRUE, max_nchar=254) {
        stopifnot(is(x, "SpatialLinesDataFrame"))
	df <- as(x, "data.frame")
	df <- data.frame(SL_ID=I(row.names(df)), df)
	pls <- .SpL2lineslist(as(x, "SpatialLines"))
	suppressWarnings(write.linelistShape(pls, df, file=fn,
	    factor2char = factor2char, max_nchar=max_nchar))
}


.shp2LinesDF <- function(shp, proj4string=CRS(as.character(NA)), IDs,
        delete_null_obj=FALSE) {
	if (class(shp) != "Map") stop("shp not a Map object")
	shp.type <- attr(shp$Shapes, "shp.type")
	if (!shp.type %in% c("arc", "poly")) 
		stop("not an arc or poly Map object")
# birds NULL part Allen H. Hurlbert 090610 copied from .Map2PolyDF
# Harlan Harris 100907
        nullParts <- sapply(shp$Shapes, function(x) x$nParts) == 0
        if (delete_null_obj) {
	    nullParts <- which(nullParts)
	    if (length(nullParts) > 0L) {
		for (i in length(nullParts):1)
	            shp$Shapes[[nullParts[i]]] <- NULL
                attr(shp$Shapes,'nshps') <- attr(shp$Shapes,'nshps') - 
                    length(nullParts)
                shp$att.data <- shp$att.data[-nullParts,]
                warning(paste("Null objects with the following", 
                    "indices deleted:", paste(nullParts, collapse=", ")))
              }
        } else {
# birds NULL part Allen H. Hurlbert 090610
# Harlan Harris 100907
            if (any(nullParts))
               stop(paste("NULL geometry found:", paste(which(nullParts),
                   collapse=", "),
                   "\n               consider using delete_null_obj=TRUE"))
	}
	df <- shp$att.data
	shapes <- shp$Shapes
	n <- length(shapes)
	LinesList <- vector(mode="list", length=n)
	if (missing(IDs)) IDs <- as.character(sapply(shapes, 
		function(x) x$shpID))
	if (length(IDs) != n) stop("IDs length differs from number of lines")
	row.names(df) <- IDs
	for (i in 1:n) {
		LinesList[[i]] <- .shapes2LinesList(shapes[[i]], ID=IDs[i])
	}
	SL <- SpatialLines(LinesList, proj4string=proj4string)
	res <- SpatialLinesDataFrame(SL, data=df)
	res
}

.shapes2LinesList <- function(shape, ID) {
	nParts <- attr(shape, "nParts")
	Pstart <- shape$Pstart
	nVerts <- nrow(shape$verts)
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
	res <- vector(mode="list", length=nParts)
	for (i in 1:nParts) {
		res[[i]] <- Line(coords=shape$verts[from[i]:to[i],,drop=FALSE])
	}
	Lines <- Lines(res, ID=ID)
	Lines
}

.SpL2lineslist <- function(x) {
	pls <- slot(x, "lines")
	n <- length(pls)
	res <- vector(mode="list", length=n)
	for (i in 1:n) {
		xyL <- lapply(slot(pls[[i]], "Lines"), 
			coordinates)
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
	}
	class(res) <- "lineslist"
	invisible(res)
}


