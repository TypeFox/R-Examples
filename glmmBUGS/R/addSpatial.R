`addSpatial` <- function(map, raggedArray=NULL, effect=NULL, prefix=NULL) {
	# find the name of the effect if it isn't specified
    if(is.null(effect)) {
    	# assume the effect is at the lowest level
    	# the lowest effect should have an Neffect and Seffect component
    	effect = grep("^(N|S)", names(raggedArray), value=T)
    	effect = substr(effect, 2,  10000)
    	effect = table(effect)
    	effect = names(effect)[effect==2]
    	if(length(effect)!= 1) {
        	warning("not sure which effect is the spatial one")
        	effect = "region"
    	}
    }
	
	# get the names of the points or regions
	if( any(names(map)==effect) ) {
		regionID = as.character(map[[effect]])
    } else 
		regionID = NULL
	
	
	spatialBugs(map, raggedArray, effect, prefix, regionID)
}   

spatialBugs = function(map, raggedArray, effect, prefix, regionID) {
	UseMethod("spatialBugs")
}


spatialBugs.SpatialPolygons <- function(map, raggedArray=NULL, effect, prefix=NULL, regionID=NULL) {
	
    map = spdep::poly2nb(map,row.names=regionID)
	spatialBugs(map, raggedArray, effect, prefix, regionID)
}

spatialBugs.nb<- function(map, raggedArray=NULL, effect, prefix=NULL, regionID=NULL) {
    if(is.null(regionID))
		regionID = attributes(map)$region.id

	map=spdep::nb2WB(map)	
    names(map$num) = regionID
	spatialBugs.list(map, raggedArray, effect, prefix, regionID)
}

spatialBugs.list = function(map, raggedArray=NULL, effect, prefix=NULL, regionID=NULL) {

  
  # make sure this is a valid winbugs adjacency list
    if(!all(c("adj", "num") %in% names(map))) {
      warning("adj or num are missing in the map")
    }
    if(is.null(regionID))
	    regionID = names(map$num)
  
  # names of the random effects
  effectNames = names(raggedArray[[paste("S", prefix, effect, sep="")]])
  effectNames = effectNames[effectNames != "end"]
  if(!length(effectNames))
    effectNames = regionID


  # check to see each effect name is a region
  if(!all(effectNames %in% regionID))
    warning("some random effects don't appear to be regions")


Sspatial = seq(1, length(regionID))
names(Sspatial) = regionID
raggedArray[[paste("Sspatial", prefix, effect, sep="")]] = Sspatial[effectNames]


  raggedArray[[paste("adj", prefix, effect,  sep="")]] = map$adj
  raggedArray[[paste("num", prefix, effect, sep="")]] = map$num
  if(is.null(map$weights)) {
    raggedArray[[paste("weights", prefix, effect, sep="")]]= rep(1, length(map$adj))
  } else {
    raggedArray[[paste("weights", prefix, effect, sep="")]] = map$weights
  }
  raggedArray[[paste("N", prefix, effect, "Spatial", sep="")]] = length(map$num)

  raggedArray

}


spatialBugs.SpatialPoints=function(map, raggedArray=NULL, effect, prefix=NULL, regionID=NULL) {
	spatialBugs.matrix(map@coords, raggedArray, effect, prefix, regionID)
}


spatialBugs.data.frame=function(map, raggedArray=NULL, effect=NULL, prefix=NULL, regionID=NULL) {
		
	if(all(c("x","y") %in% names(map))) {
		map = map[,c("x","y")]
	} else if(dim(map)[2] == 2) {
		names(map) = c("x","y")
	} else {
		warning("can't find coordinate columns, should be named x and y")
		
	}

	spatialBugs.matrix(as.matrix(map), raggedArray, effect, prefix, regionID)
	
}

spatialBugs.complex=function(map, raggedArray=NULL, effect=NULL, prefix=NULL, regionID=NULL) {
	if(is.null(regionID))
		regionID = names(map)
	
	spatialBugs(cbind(x=Re(map), y=Im(map)), raggedArray, effect, prefix, regionID)
}



spatialBugs.matrix=function(map, raggedArray=NULL, effect=NULL, prefix=NULL, regionID=NULL) {

	if(all(c("x","y") %in% colnames(map))) {
		map = map[,c("x","y")]
	} else if(dim(map)[2] == 2) {
		colnames(map) = c("x","y")
	} else {
		warning("can't find coordinate columns, should be named x and y")
	}
	
	
	effect = paste(prefix, effect, sep="")

	# reorder the rows
	if(!is.null(regionID))
		rownames(map) = regionID
	theNames = names(raggedArray[[paste("S", effect, sep="")]])
	theNames = theNames[theNames != "end"]
	if(!all(theNames %in% rownames(map)))
		warning("can't find some locations with data in the coordinates specified")
		
	if(!is.null(rownames(map))) {
		map = map[theNames,]
	}
	
	# add to ragged array
	raggedArray[[paste("xSpatial", effect, sep="")]] = map[,"x"]
	raggedArray[[paste("ySpatial", effect, sep="")]] = map[,"y"]
	
	raggedArray[[paste("N", effect, "Spatial", sep="")]] = dim(map)[1]

	raggedArray
}
