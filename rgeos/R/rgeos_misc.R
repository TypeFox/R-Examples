RGEOSMiscFunc = function(spgeom, byid, func) {
    stopifnot(is.logical(byid))
    if (!is.na(is.projected(spgeom)) && !is.projected(spgeom))
     warning("Spatial object is not projected; GEOS expects planar coordinates")
    byid = as.logical(byid)
    if (is.na(byid)) stop("Invalid value for byid, must be logical")
    if (inherits(spgeom, "SpatialPolygons") && get_do_poly_check() && notAllComments(spgeom)) 
        spgeom <- createSPComment(spgeom)
    if (func == "rgeos_area")    
        x <- .Call("rgeos_area", .RGEOS_HANDLE, spgeom, byid, PACKAGE="rgeos")
    else if (func == "rgeos_length")    
        x <- .Call("rgeos_length", .RGEOS_HANDLE, spgeom, byid, PACKAGE="rgeos")
    else stop("no such function:", func)

    if(byid) names(x) <- unique(row.names(spgeom))
    
    return(x)
}

gArea = function(spgeom, byid=FALSE) {
    return( RGEOSMiscFunc(spgeom,byid,"rgeos_area") )
}

gLength = function(spgeom, byid=FALSE) {
    return(RGEOSMiscFunc(spgeom,byid,"rgeos_length"))
}


RGEOSDistanceFunc = function(spgeom1, spgeom2, byid, func, densifyFrac = 1) {
    stopifnot(is.logical(byid))
    if (!is.na(is.projected(spgeom1)) && !is.projected(spgeom1))
     warning("Spatial object 1 is not projected; GEOS expects planar coordinates")
    if (!is.null(spgeom2) && !is.na(is.projected(spgeom2)) && 
        !is.projected(spgeom2))
     warning("Spatial object 2 is not projected; GEOS expects planar coordinates")
    byid = as.logical(byid)
    if (any(is.na(byid)) ) stop("Invalid value for byid, must be logical")

    if( length(byid) < 1 || length(byid) > 2 )
        stop("Invalid length for byid, must be of length 1 or 2")

    if (length(byid) == 1)
        byid <- rep(byid,2)

    if(!is.null(spgeom1) & !is.null(spgeom2)) {
        if(!identical(spgeom1@proj4string,spgeom2@proj4string))
            warning("spgeom1 and spgeom2 have different proj4 strings")
    }

    if (func == "rgeos_hausdorffdistancedensify")
        x <- .Call("rgeos_hausdorffdistancedensify", .RGEOS_HANDLE, spgeom1, spgeom2, densifyFrac, byid, PACKAGE="rgeos")
    else if (func == "rgeos_hausdorffdistance")
        x <- .Call("rgeos_hausdorffdistance", .RGEOS_HANDLE, spgeom1, spgeom2, byid, PACKAGE="rgeos")
    else if (func == "rgeos_distance")
        x <- .Call("rgeos_distance", .RGEOS_HANDLE, spgeom1, spgeom2, byid, PACKAGE="rgeos")
    else stop("no such function:", func)

    
    if(any(byid)) {
        id1 = row.names(spgeom1) 
        if (is.null(spgeom2))
            id2 = id1
        else 
            id2 = row.names(spgeom2)
        
        if (length(id1) == ncol(x)) colnames(x) = id1
        if (length(id2) == nrow(x)) rownames(x) = id2
    }

    return(x)
} 

gDistance = function(spgeom1, spgeom2=NULL, byid=FALSE, hausdorff=FALSE, densifyFrac = NULL) {
	if (hausdorff) {
		if(is.null(densifyFrac)) {
			return( RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_hausdorffdistance") )
		} else {
			if (!is.numeric(densifyFrac))
		        stop("densifyFrac must be numeric")

		    if (densifyFrac > 1 | densifyFrac <= 0)
		        stop("densifyFrac must be in the range (0,1]")

		    return( RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_hausdorffdistancedensify", densifyFrac) )
		}
	} else {
    	return( RGEOSDistanceFunc(spgeom1, spgeom2, byid, "rgeos_distance") )
	}
}

gWithinDistance = function(spgeom1, spgeom2=NULL, dist, byid=FALSE, hausdorff=FALSE, densifyFrac=NULL) {
    #TODO - include a tolerance?
    return( gDistance(spgeom1,spgeom2,byid, hausdorff, densifyFrac) <= dist )
}


gTopoDim <- function(obj) {
    td <- NULL
    if (inherits(obj, "SpatialPolygons")) td <- 2L
    if (inherits(obj, "SpatialRings")) td <- 2L
    if (inherits(obj, "SpatialLines")) td <- 1L
    if (inherits(obj, "SpatialPoints")) td <- 0L
    if (is.null(td)) stop("class not supported:", class(obj))
    td
}

#Deprecated function names
RGEOSArea = function(spgeom, byid=FALSE) {
    .Deprecated("gArea")
    return( gArea(spgeom, byid) )
}
RGEOSLength = function(spgeom, byid=FALSE) {
    .Deprecated("gLength")
    return( gLength(spgeom, byid) )
}
RGEOSDistance = function(spgeom1, spgeom2=NULL, byid=FALSE) {
    .Deprecated("gDistance")
    return( gDistance(spgeom1, spgeom2, byid) )
}
RGEOSisWithinDistance = function(spgeom1, spgeom2=NULL, dist, byid=FALSE) {
    .Deprecated("gWithinDistance")
    return( gWithinDistance(spgeom1, spgeom2, dist, byid) )
}
RGEOSHausdorffDistance = function(spgeom1, spgeom2=NULL, byid=FALSE) {
    .Deprecated("gDistance")
    return( gDistance(spgeom1, spgeom2, byid, TRUE) )
}
