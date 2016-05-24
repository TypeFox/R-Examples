#Date: 20.07.2013

# Identying the closer stretch(es) to a station	
# For SpPolygons and SpPoints intersection just work for length(SpPoints) = 1
	
RiverStation = function(x, y, window=100){
	# x = SpatialPointsDataFrame, length =1
	# y = SpatialLinesDataFrame
	# window = size of the square (window) around the point
	# RiverStation 
	
	# creating a square around the station
		# x = station
		# y = riverIO
		cx = slot(x, "coords")[1]
		cy = slot(x, "coords")[2]
		d = window/2 
		sr1 = Polygon(cbind(c(cx-d, cx+d, cx+d, cx-d, cx-d), c(cy-d, cy-d, cy+d, cy+d, cy-d)))
		srs1 = Polygons(list(sr1), "s1")
		sp = SpatialPolygons(list(srs1))
     slot(sp, "proj4string") = slot(y, "proj4string")
			#plot(sp)
			#plot(sp, add = T)
			#plot(station, add=T)
			   
	# intersection of the square and the rivers
		#str(river, level=1)
   	#str(sp)
   	idSq = gIntersects(y, sp, byid=T); idSq
   	id = 
   	if(length(idSq)==length(which(idSq==FALSE))){ # there are not intersection
   		return(riverStation=0)
   	}else{
	   	#class(idSq)
   		idSq = which(idSq == "TRUE"); idSq
   		#class(idSq)
   		riverStation = y[idSq,]
   		return(riverStation)
   	}
}