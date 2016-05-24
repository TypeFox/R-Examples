sp2tmap <- function(SP) {
	if (!inherits(SP, "SpatialPolygons"))
		stop("not a SpatialPolygons object")
	pls <- slot(SP, "polygons")
	IDs <- sapply(pls, function(x) slot(x, "ID"))
	n <- length(IDs)
	iIDs <- as.integer(1:n)
	cID <- NULL
	cX <- NULL
	cY <- NULL
	for (i in iIDs) {
		pl <- slot(pls[[i]], "Polygons")
		m <- length(pl)
		for (j in 1:m) {
			crds <- slot(pl[[j]], "coords")
			if (is.null(cID)) { 
				cID <- i
				cX <- as.numeric(NA)
				cY=as.numeric(NA)
			} else { 
				cID <- c(cID, i)
				cX <- c(cX, as.numeric(NA))
				cY <- c(cY, as.numeric(NA))
			}
			cID <- c(cID, rep(i, nrow(crds)))
			cX <- c(cX, crds[,1])
			cY <- c(cY, crds[,2])
		}
	}
	res <- data.frame("_ID"=cID, "_X"=cX, "_Y"=cY, check.names=FALSE)
	names(iIDs) <- IDs
	attr(res, "ID_names") <- iIDs
	res
}
