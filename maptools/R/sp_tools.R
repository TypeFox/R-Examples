MapGen2SL <- function(file, proj4string=CRS(as.character(NA))) {
	con <- file(file, "r")
	hold <- readLines(con)
	close(con)
	if (length(hold) == 500000L) warning("500,000 point limit reached")
	starts <- which(hold == "# -b")
	n <- length(starts)
	if (n < 1) stop("Not a Mapgen format file")
	res <- vector(mode="list", length=n)
	IDs <- paste("L", 1:n, sep="_")
	for (i in 1:n) {
		if (i < n) {
			x <- t(sapply(strsplit(hold[(starts[i]+1):
				(starts[i+1]-1)], "\t"), as.numeric))
		} else {
			x <- t(sapply(strsplit(hold[(starts[i]+1):
				length(hold)], "\t"), as.numeric))
		}
		res[[i]] <- Lines(list(Line(x
#, proj4string=proj4string
)), ID=IDs[i])
	}
	SL <- SpatialLines(res, proj4string=proj4string)
	SL
}

ArcObj2SLDF <- function(arc, proj4string=CRS(as.character(NA)), IDs) {
	df <- data.frame(arc[[1]])
	n <- length(arc[[2]])
	LinesList <- vector(mode="list", length=n)
	if (missing(IDs)) IDs <- paste("L", 1:n, sep="_")
	if (length(IDs) != n) stop("IDs length differs from number of arcs")
	row.names(df) <- IDs
	for (i in 1:n) {
		crds <- cbind(arc[[2]][[i]][[1]], arc[[2]][[i]][[2]])
		LinesList[[i]] <- Lines(list(Line(coords=crds
#, proj4string=proj4string
)), ID=IDs[i])
	}
	SL <- SpatialLines(LinesList, proj4string=proj4string)
	res <- SpatialLinesDataFrame(SL, data=df)
	res
}

ContourLines2SLDF <- function(cL, proj4string=CRS(as.character(NA))) {
	if (length(cL) < 1L) stop("cL too short")
	cLstack <- tapply(1:length(cL), sapply(cL, function(x) x[[1]]), 
		function(x) x, simplify=FALSE)
	df <- data.frame(level=names(cLstack))
	m <- length(cLstack)
	res <- vector(mode="list", length=m)
	IDs <- paste("C", 1:m, sep="_")
	row.names(df) <- IDs
	for (i in 1:m) {
		res[[i]] <- Lines(.contourLines2LineList(cL[cLstack[[i]]]#, 
#			proj4string=proj4string
), ID=IDs[i])
	}
	SL <- SpatialLines(res, proj4string=proj4string)
	res <- SpatialLinesDataFrame(SL, data=df)
	res
}
.contourLines2LineList <- function(cL#, proj4string=CRS(as.character(NA))
) {
	n <- length(cL)
	res <- vector(mode="list", length=n)
	for (i in 1:n) {
		crds <- cbind(cL[[i]][[2]], cL[[i]][[3]])
		res[[i]] <- Line(coords=crds#, proj4string=proj4string
)
	}
	res
}

# to be moved to glue with RarcInfo:

pal2SpatialPolygons <- function(arc, pal, IDs, dropPoly1=TRUE, 
	proj4string=CRS(as.character(NA))) {
	if (missing(IDs)) stop("IDs required")
	if (dropPoly1) pale <- lapply(pal[[2]][-1], function(x) x[[1]])
	else pale <- lapply(pal[[2]], function(x) x[[1]])
	if (length(pale) != length(IDs)) stop("map and IDs differ in length")
	tab <- table(factor(IDs))
	n <- length(tab)
	IDss <- names(tab)
	reg <- match(IDs, IDss)
	belongs <- lapply(1:n, function(x) which(x == reg))
# assemble the list of Polygons
	Srl <- vector(mode="list", length=n)
	for (i in 1:n) {
		bi <- belongs[[i]]
		nParts <- length(bi)
		palei_list <- list()
		for (j in 1:nParts) {
			this <- bi[j]
			paleij <- pale[[this]]
			if (any(paleij == 0)) {
				zeros <- which(paleij == 0)
				palei_list <- c(palei_list, 
					list(paleij[1:(zeros[1]-1)]))
				for (k in 1:length(zeros)) {
					if (k == length(zeros)) {
						lp <- length(paleij)
						lz <- zeros[length(zeros)]
						palei_list <- c(palei_list, 
						    list(paleij[(lz+1):lp]))
					} else {
						zk <- zeros[k]
						zk1 <- zeros[k+1]
						palei_list <- c(palei_list, 
						    list(paleij[(zk+1):(zk1-1)]))
					}
				}
			} else palei_list <- c(palei_list, list(paleij))
		}
		nParts <- length(palei_list)
		srl <- vector(mode="list", length=nParts)
		for (j in 1:nParts) {
			paleij <- palei_list[[j]]
			nArcs <- length(paleij)
			x <- NULL
			y <- NULL
			for (k in 1:nArcs) {
				kk <- paleij[k]
				if (kk > 0) {
					x <- c(x, arc[[2]][[kk]][[1]])
					y <- c(y, arc[[2]][[kk]][[2]])
				} else {
					x <- c(x, rev(arc[[2]][[-kk]][[1]]))
					y <- c(y, rev(arc[[2]][[-kk]][[2]]))
				}
			}
			if ((x[1] != x[length(x)]) || (y[1] != y[length(y)])) {
				x <- c(x, x[1])
				y <- c(y, y[1])
			}
			srl[[j]] <- Polygon(coords=cbind(x, y))	
		}
		Srl[[i]] <- Polygons(srl, ID=IDss[i])
	}
	res <- as.SpatialPolygons.PolygonsList(Srl, proj4string=proj4string)
	res
}


