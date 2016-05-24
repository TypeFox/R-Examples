gSimplify = function(spgeom, tol, topologyPreserve=FALSE) {

	getCutEdges = as.logical(topologyPreserve)
	if (is.na(topologyPreserve))
		stop("Invalid value for topologyPreserve, must be logical")
	
    if (inherits(spgeom, "SpatialPolygons") && get_do_poly_check() && notAllComments(spgeom)) 
        spgeom <- createSPComment(spgeom)
    id = row.names(spgeom)
    return( .Call("rgeos_simplify", .RGEOS_HANDLE, spgeom, tol, id, 
									FALSE, topologyPreserve, PACKAGE="rgeos") )
}

gPolygonize = function( splist, getCutEdges=FALSE) {
	
	if (!is.list(splist))
		splist = list(splist)
        if (!all(sapply(splist, inherits, "SpatialLines")))
            stop("list of SpatialLines object required")

	p4slist = lapply(splist,function(x) x@proj4string)
#        splist <- lapply(splist, function(s) {
#            if (inherits(s, "SpatialPolygons") && get_do_poly_check() && notAllComments(s)) {
#                createSPComment(s)
#            } else {
#                s
#            }
#        })
	
	p4s = p4slist[[1]]
	if (length(p4slist) != 1) {
		for(i in 2:length(p4slist)) {
			if (!identical(p4s, p4slist[[i]]))
				stop("spgeoms in splist have different proj4strings")
		}
	}
	

	getCutEdges = as.logical(getCutEdges)
	if (is.na(getCutEdges))
		stop("Invalid value for getCutEdges, must be logical")
	
	nid = sum(sapply(splist, function(x) sapply(slot(x, "lines"),
            function(y) length(slot(y, "Lines")))))
#           sum(sapply(splist,function(x) length(unlist(row.names(x)))))
    id = as.character(1:nid)

    return( .Call("rgeos_polygonize", .RGEOS_HANDLE, splist, id, p4s,
									getCutEdges, PACKAGE="rgeos") )
}





TopologyFunc = function(spgeom, id, byid, func) {
    
    stopifnot(is.logical(byid))
    byid = as.logical(byid)
    if (is.na(byid)) 
        stop("Invalid value for byid, must be logical")
    
    curids = unique(row.names(spgeom))
    if (is.null(id)) {
        if (byid)   id = curids
        else        id = "1"
    }
    id = as.character(id)

    if (inherits(spgeom, "SpatialPolygons") && get_do_poly_check() && notAllComments(spgeom)) spgeom <- createSPComment(spgeom)
    
    if ( length(id) != length(unique(id)) )
        stop("Non-unique values for id ")
    
    if ( !(!byid && length(id) == 1) && !(byid && length(id) == length(curids)) )
        stop("Invalid number of values in id" ) 
    if (func == "rgeos_envelope")
        x <- .Call("rgeos_envelope", .RGEOS_HANDLE, spgeom, id, byid, PACKAGE="rgeos")
    else if (func == "rgeos_convexhull")
        x <- .Call("rgeos_convexhull", .RGEOS_HANDLE, spgeom, id, byid, PACKAGE="rgeos")
    else if (func == "rgeos_boundary")
        x <- .Call("rgeos_boundary", .RGEOS_HANDLE, spgeom, id, byid, PACKAGE="rgeos")
    else if (func == "rgeos_getcentroid")
        x <- .Call("rgeos_getcentroid", .RGEOS_HANDLE, spgeom, id, byid, PACKAGE="rgeos")
    else if (func == "rgeos_pointonsurface")
        x <- .Call("rgeos_pointonsurface", .RGEOS_HANDLE, spgeom, id, byid, PACKAGE="rgeos")
    else if (func == "rgeos_unioncascaded")
        x <- .Call("rgeos_unioncascaded", .RGEOS_HANDLE, spgeom, id, byid, PACKAGE="rgeos")
    else stop("no such function:", func)


    return( x )
}

gEnvelope = function(spgeom, byid=FALSE, id = NULL) {
    return( TopologyFunc(spgeom,id,byid,"rgeos_envelope") ) 
}
gConvexHull = function(spgeom, byid=FALSE, id = NULL) {
    return( TopologyFunc(spgeom,id,byid,"rgeos_convexhull") ) 
}
gBoundary = function(spgeom, byid=FALSE, id = NULL) {
     return( TopologyFunc(spgeom,id,byid,"rgeos_boundary") ) 
}
gCentroid = function(spgeom, byid=FALSE, id = NULL) {
    return( TopologyFunc(spgeom,id,byid,"rgeos_getcentroid") ) 
}
gPointOnSurface = function(spgeom, byid=FALSE, id = NULL) {
    return( TopologyFunc(spgeom,id,byid,"rgeos_pointonsurface") ) 
}
gLineMerge = function(spgeom, byid=FALSE, id = NULL) {
#    return( TopologyFunc(spgeom,id,byid,"rgeos_linemerge") ) 
    if (!inherits(spgeom,"SpatialLines"))
        stop("Invalid geometry, may only be applied to lines")
    spgeom <- as(spgeom, "SpatialLines")
    if (is.null(id))
        id = rep("1",length(row.names(spgeom)))

#    if (any(is.na(id))) stop("No NAs permitted in id")

    ids <- split(1:length(id), id)
    out <- vector(mode="list", length=length(ids))
    for (i in seq(along=ids)) {
        out[[i]] <- .Call("rgeos_linemerge", .RGEOS_HANDLE,
        spgeom[ids[[i]]], names(ids)[i], FALSE, PACKAGE="rgeos") 
    }
    res <- do.call("rbind.SpatialLines", out)
    res
}

gUnionCascaded = function(spgeom, id = NULL) {
    
    if (!inherits(spgeom,"SpatialPolygons"))
        stop("Invalid geometry, may only be applied to polygons")
    spgeom <- as(spgeom, "SpatialPolygons")

    if (is.null(id))
        id = rep("1",length(row.names(spgeom)))

#    if (any(is.na(id))) stop("No NAs permitted in id")

    if (get_do_poly_check() && notAllComments(spgeom)) spgeom <- createSPComment(spgeom)

    ids <- split(1:length(id), id)
    sl <- sapply(ids, length)
    out <- vector(mode="list", length=length(ids))
    for (i in seq(along=ids)) {
        out[[i]] <- TopologyFunc(groupID(spgeom[ids[[i]]], id[ids[[i]]]),
            names(ids)[i], TRUE, "rgeos_unioncascaded")
    }
    res <- do.call("rbind.SpatialPolygons", out)

    res
}

gUnaryUnion = function(spgeom, id = NULL) {

    if (version_GEOS0() < "3.3.0")
        stop("No UnaryUnion in this version of GEOS")
    
    if (!inherits(spgeom,"SpatialPolygons"))
        stop("Invalid geometry, may only be applied to polygons")
    spgeom <- as(spgeom, "SpatialPolygons")
    if (is.null(id))
        id = rep("1",length(row.names(spgeom)))

#    if (any(is.na(id))) stop("No NAs permitted in id")

    if (get_do_poly_check() && notAllComments(spgeom)) spgeom <- createSPComment(spgeom)

    ids <- split(1:length(id), id)
    out <- vector(mode="list", length=length(ids))
    for (i in seq(along=ids)) {
        out[[i]] <- .Call("rgeos_unaryunion", .RGEOS_HANDLE,
        spgeom[ids[[i]]], names(ids)[i], FALSE, PACKAGE="rgeos") 
    }
    res <- do.call("rbind.SpatialPolygons", out)
    res
}

gDelaunayTriangulation <- function(spgeom, tolerance=0.0, onlyEdges=FALSE) {

    if (version_GEOS0() < "3.4.0")
        stop("No DelaunayTriangulation in this version of GEOS")

    if (!inherits(spgeom, "SpatialPoints"))
        stop("Invalid geometry, may only be applied to points")
    if (nrow(zerodist(spgeom)) > 0)
        stop("duplicate points not permitted")
    stopifnot(is.numeric(tolerance))
    stopifnot(length(tolerance) == 1)
    stopifnot(is.logical(onlyEdges))
    stopifnot(length(onlyEdges) == 1)

    .Call("rgeos_delaunaytriangulation", .RGEOS_HANDLE,
        spgeom, tolerance, onlyEdges, PACKAGE="rgeos")
    
}

gNode <- function(spgeom){

    if (version_GEOS0() < "3.4.0")
        stop("No noding in this version of GEOS")

    if (!inherits(spgeom, "SpatialLines"))
        stop("Invalid geometry, may only be applied to lines")

    .Call("rgeos_node", .RGEOS_HANDLE, spgeom, PACKAGE="rgeos")

}

RGEOSEnvelope = function(spgeom, byid=FALSE, id = NULL) {
    .Deprecated("gEnvelope")
    return( gEnvelope(spgeom, id, byid) )
}
RGEOSConvexHull = function(spgeom, byid=FALSE, id = NULL) {
    .Deprecated("gConvexHull")
    return( gConvexHull(spgeom, id, byid) )
}
RGEOSBoundary = function(spgeom, byid=FALSE, id = NULL) {
    .Deprecated("gBoundary")
    return( gBoundary(spgeom, id, byid) )
}
RGEOSGetCentroid = function(spgeom, byid=FALSE, id = NULL) {
    .Deprecated("gCentroid")
    return( gCentroid(spgeom, id, byid) )
}
RGEOSPointOnSurface = function(spgeom, byid=FALSE, id = NULL) {
    .Deprecated("gPointOnSurface")
    return( gPointOnSurface(spgeom, id, byid) )
}
RGEOSLineMerge = function(spgeom, byid=FALSE, id = NULL) {
    .Deprecated("gLineMerge")
    return( gLineMerge(spgeom, id, byid) )
}
RGEOSUnionCascaded = function(spgeom, id = NULL) {
    .Deprecated("gUnionCascaded")
    return( gUnionCascaded(spgeom, id) )
}




