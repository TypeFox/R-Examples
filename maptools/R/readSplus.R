#
#Read maps in S-plus format (as exported by WinBUGS)
#
readSplus<-function(file, proj4string=CRS(as.character(NA)) )
{

	lns<-readLines(file)
	nlines<-length(lns)
	nareas<-as.integer(strsplit(lns[1], ":")[[1]][2])

	offset<-1
	if(lns[2]==""){offset<-offset+1}

        xscale <- 1
        yscale <- 1
        scales <- grep("cale", lns[(offset+1):(offset+2)])
	if (length(scales) > 0L) {
            if (length(scales) < 2L) 
                stop("Only one scale given")
            xsc <- grep("x|X", lns[(offset+1):(offset+2)])
            ysc <- grep("y|Y", lns[(offset+1):(offset+2)])
            xscale <- strsplit(lns[(offset+xsc)], ":| |\t|\r\n")[[1]]
            if (any(xscale == "")) xscale <- xscale[-which(xscale == "")]
            xscale <- as.numeric(xscale[2])
            yscale <- strsplit(lns[(offset+xsc)], ":| |\t|\r\n")[[1]]
            if (any(yscale == "")) yscale <- yscale[-which(yscale == "")]
            yscale <- as.numeric(yscale[2])
            offset <- offset+2
	    if(lns[(offset+1)]==""){offset<-offset+1}
        }

	IDs<-lapply(lns[offset+1:nareas], function(X){strsplit(X, " |\t|\r\n")})
	IDs<-matrix(unlist(IDs), ncol=2, byrow=TRUE)

	offset<-offset+nareas

	if(lns[offset+1]==""){offset<-offset+1}

        END <- which(lns == "END")
        lns <- lns[(offset+1):(END-1)]
        if (any(lns == "")) {
            empty <- which(lns == "")
            lns <- lns[-empty]
        }

#	polys<-read.table(file, skip=offset, nrows=nlines-offset-END)
        polys <- lapply(lns, function(X){strsplit(X, " |\t|\r\n")})
        polys <- matrix(unlist(polys), ncol=3, byrow=TRUE)
        wNA <- which(polys[,1] == "NA")
        is.na(polys[wNA,1]) <- TRUE
        is.na(polys[wNA,2]) <- TRUE
        is.na(polys[wNA,3]) <- TRUE

	polys2<-cbind(xscale*as.numeric(polys[,2]),
            yscale*as.numeric(polys[,3]))

	lpolys<-.NAmat2xyList(polys2)
#	llpolys<-unlist(lapply(lpolys, nrow))

	#OpenBUGS seems to put a line with NAs just before the END
	#We need to remove it
	if(is.na(polys[nrow(polys),1]))
		polys<-polys[-nrow(polys),]
		

#	idx<- c(1, cumsum(2+llpolys[-length(lpolys)]))
        wNA <- which(is.na(polys[, 1]))


        idx <- c(1, wNA+1)
	polysIDs<-polys[idx, 1]

	belongs<-lapply(1:nareas, function(i){which(polysIDs==IDs[i,2])})

	Srl <- vector(mode = "list", length = nareas)

	for (i in 1:nareas) {
		nParts <- length(belongs[[i]])
		srl <- vector(mode = "list", length = nParts)
		for (j in 1:nParts) {
                        crds <- lpolys[[belongs[[i]][j]]]
                        nc <- nrow(crds)
                        if (crds[1,1] != crds[nc,1] || crds[1,2] != crds[nc,2])
                            crds <- rbind(crds, crds[1,,drop=FALSE])
			srl[[j]] <- Polygon(coords = crds)
		}
		Srl[[i]] <- Polygons(srl, ID = IDs[i,2])
	}

	res <- as.SpatialPolygons.PolygonsList(Srl, proj4string = proj4string)
	res
}

