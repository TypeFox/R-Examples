###Redifining spTransform, because it changes the class of the object to SpatialPointsDataFrame 
setMethod(f = "spTransform", 
	  signature = c(x = ".MoveTrack", CRSobj = "missing"), 
	  function(x, center=FALSE, ...){
		  if(!center)
			  stop('spTransform without center or proj string does not make much sense')
		  spTransform(x=x, center=center, CRSobj="+proj=aeqd +ellps=WGS84")
	  })

setMethod(f = "spTransform", 
	  signature = c(x = ".MoveTrack", CRSobj = "character"), 
	  function(x, CRSobj, center=FALSE, ...){
		  spTransform(x=x, CRSobj=CRS(CRSobj), center=center)
	  })

setMethod(f = "spTransform", 
	  signature = c(x = ".MoveTrack", CRSobj = "CRS"), 
	  function(x, CRSobj, center=FALSE, ...){
		  if (center){ 
			  if(!isLonLat(x)){
				  stop('Center only works with Longitude Latitude projection')
			  }
			  mid.range.lon <- (max(coordinates(x)[ ,1])+min(coordinates(x)[ ,1]))/2
			  mid.range.lat  <- (max(coordinates(x)[ ,2])+min(coordinates(x)[ ,2]))/2
			  CRSobj <- CRS(paste(CRSobj@projargs," +lon_0=",mid.range.lon," +lat_0=", mid.range.lat, sep=""))
		  } 

		  coordsnew <- spTransform(x=as(x,'SpatialPointsDataFrame'), CRSobj=CRSobj)
		  x <- new(class(x), coordsnew,x )
		  return(x)
	  })
