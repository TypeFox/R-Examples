trans.mat <- function(bathy,min.depth=0,max.depth=NULL) {
	
	# require(gdistance)
	
	ras <- bathy
	ras[bathy > min.depth] <- 0.00000001
	ras[bathy <= min.depth] <- 1
	if (!is.null(max.depth)) ras[bathy <= max.depth] <- 0.00000001
	
	lat <- as.numeric(colnames(bathy))
	lon <- as.numeric(rownames(bathy))
	
	r <- raster::raster(ncol=nrow(bathy),nrow=ncol(bathy),xmn=min(lon),xmx=max(lon),ymn=min(lat),ymx=max(lat))
	raster::values(r) <- as.vector(ras[,rev(1:ncol(ras))])
	
	trans <- gdistance::transition(r, transitionFunction = mean, directions = 16)
	transC <- gdistance::geoCorrection(trans,type="c",multpl=FALSE)

	return(transC)
}