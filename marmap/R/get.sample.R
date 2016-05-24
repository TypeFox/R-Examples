get.sample=function(mat, sample, col.lon, col.lat, ...){
	
	locator(n=2,type="n")->coord
	as.numeric(rownames(mat), na.rm=TRUE) -> lon
	as.numeric(colnames(mat), na.rm=TRUE) -> lat
	
	if(length(coord$x) == 1) {	
		warning("Please choose two points from the map")
		}
		
	if(length(coord$x) == 2) {
		
		rect(min(coord$x),min(coord$y),max(coord$x),max(coord$y),...)
		
		which(abs(lon-coord$x[1])==min(abs(lon-coord$x[1]))) -> x1
		which(abs(lat-coord$y[1])==min(abs(lat-coord$y[1]))) -> y1
		which(abs(lon-coord$x[2])==min(abs(lon-coord$x[2]))) -> x2
		which(abs(lat-coord$y[2])==min(abs(lat-coord$y[2]))) -> y2

		new.bathy = mat[x1:x2, y1:y2]		
		as.numeric(rownames(new.bathy), na.rm=TRUE) -> lon2
		as.numeric(colnames(new.bathy), na.rm=TRUE) -> lat2
		
		range(lon2) -> rlon2
		range(lat2) -> rlat2			
		
		subset(sample, sample[,col.lon]>rlon2[1] & sample[,col.lat]>rlat2[1]) -> s1
		subset(s1, s1[,col.lon]<rlon2[2] & s1[,col.lat]<rlat2[2]) -> s2
		
			if(length(s2[,1]) == 0) cat("No sample in the selected region")
			if(length(s2[,1]) != 0) return(droplevels(s2))
		
		
		}

}