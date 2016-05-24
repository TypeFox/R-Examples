trimdata <-
function(rdata, map, Xcol=2, Ycol=3, rectangle=F, buffer=0.05) {
	if(class(map)=="SpatialPolygonsDataFrame") {
		mappoly=SpatialPolygons2PolySet(map) 
		mr = c(range(mappoly$X),range(mappoly$Y))	# map range
	} else 
	if (class(map)=="map") {
		map_sp = map2SpatialLines(map,proj4string=CRS("+proj=longlat +datum=WGS84"))	
		mappoly=SpatialLines2PolySet(map_sp)
		mr = map$range		# map range
	} else
		stop(paste("map class not recognized by trimdata function"))
	# if rdata has only 2 columns assign them as X and Y
	if (length(rdata)==2) {
		Xcol = 1
		Ycol = 2
	}
    dataXmin=min(rdata[,Xcol])
	dataXmax=max(rdata[,Xcol])
	dataYmin=min(rdata[,Ycol])
	dataYmax=max(rdata[,Ycol])
	dr = c(dataXmin,dataXmax,dataYmin,dataYmax)		# data range
	
	# Print warning if the map centroid is not within the data range
	# only works if PBSmapping library is loaded
	if (requireNamespace("PBSmapping", quietly = TRUE)) {
      centroid = PBSmapping::calcCentroid(PBSmapping::combinePolys(mappoly),rollup=1)[-1]
	  mapmean=(dataXmin<centroid$X)+(centroid$X<dataXmax)+(dataYmin<centroid$Y)+(centroid$Y<dataYmax)
	  if (mapmean!=4) {				
		warning("The map centroid is not in the range of the data.  This might indicate differing projections.")
		cat(paste(c("Data X: "," to ","; Data Y: "," to "),signif(dr,digits=4),sep="",collapse=""),fill=TRUE)
		cat(paste(c("Map X: "," to ","; Map Y: "," to "),signif(mr,digits=4),sep="",collapse=""),fill=TRUE)
		cat(paste(c("Map centroid: ",", "),signif(centroid,digits=4),sep="",collapse=""),fill=TRUE)
	  }
	}
	if (rectangle==T) { 
		# Keep subset of data within a rectangular boundary using map coordinate ranges
		if (length(buffer)<=2) buffer = rep(c(mr[2]-mr[1],mr[4]-mr[3])*buffer,each=2) 
		mr = mr + c(-1,1,-1,1)*buffer
		rdata = rdata[rdata[,2] > mr[1] & rdata[,2] < mr[2] & rdata[,3] > mr[3] & rdata[,3] < mr[4],]
		cat(paste(c("Data trimmed to rectangle on X: "," to "," and Y: "," to "),
				signif(mr,digits=4),sep="",collapse=""),fill=TRUE)
	} else {
		# Keep subset of data within strict boundaries of map
		rdn = names(rdata)
		rdata$EID = 1:length(rdata[,1])
		names(rdata)[c(Xcol,Ycol)] = c("X","Y")
		maped=PBSmapping::as.EventData(rdata)
		mapfind=PBSmapping::findPolys(maped,mappoly)$EID	
		rdata = rdata[rdata$EID %in% mapfind,]
		rdata = rdata[,!(names(rdata) %in% "EID")]	# keep all but EID variable
		names(rdata) = rdn							# restore original variable names
	}
	return(rdata)
}
