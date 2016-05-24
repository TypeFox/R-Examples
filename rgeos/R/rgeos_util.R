poly_findInBoxGEOS <- function(spl, as_points=TRUE) {
    stopifnot(is(spl, "SpatialPolygons"))
    stopifnot(is.logical(as_points))
    stopifnot(!is.na(as_points))
    pls <- slot(spl, "polygons")
    res <- .Call("rgeos_poly_findInBox", .RGEOS_HANDLE, pls, as_points,
       PACKAGE="rgeos")
    res
}

gUnarySTRtreeQuery <- function(obj) {
    if(inherits(obj, "SpatialLines")) type <- "line"
    else if(inherits(obj, "SpatialPolygons")) type <- "poly"
    else if(inherits(obj, "Polygons")) type <- "Poly"
    else stop(paste("unsupported class:", class(obj)))
    if (type == "line") lst <- slot(obj, "lines")
    else if (type == "poly") lst <- slot(obj, "polygons")
    else lst <- slot(obj, "Polygons")
    res <- .Call("rgeos_unary_STRtree_query", .RGEOS_HANDLE, lst,
        PACKAGE="rgeos")
    res
}

gBinarySTRtreeQuery <- function(obj1, obj2) {
    if(inherits(obj1, "SpatialLines")) type1 <- "line"
    else if(inherits(obj1, "SpatialPolygons")) type1 <- "poly"
    else stop(paste("unsupported class:", class(obj1)))
    if(inherits(obj2, "SpatialLines")) type2 <- "line"
    else if(inherits(obj2, "SpatialPoints")) type2 <- "pts"
    else if(inherits(obj2, "SpatialPolygons")) type2 <- "poly"
    else stop(paste("unsupported class:", class(obj2)))
    if (type1 == "line") lst1 <- slot(obj1, "lines")
    else lst1 <- slot(obj1, "polygons")
    if (type2 == "line") lst2 <- slot(obj2, "lines")
    else if (type2 == "pts") lst2 <- as(obj2, "SpatialPoints")
    else lst2 <- slot(obj2, "polygons")
    res <- .Call("rgeos_binary_STRtree_query", .RGEOS_HANDLE, lst1, lst2,
        PACKAGE="rgeos")
    res
}


set_do_poly_check <- function(value) {
  stopifnot(is.logical(value))
  stopifnot(length(value) == 1)
  assign("do_poly_check", value, envir=.RGEOS_HANDLE)
}

get_do_poly_check <- function() {
  get("do_poly_check", envir=.RGEOS_HANDLE)
}

notAllComments <- function(spgeom) {
    if (!get_do_poly_check()) return(FALSE)
    if (!inherits(spgeom, "SpatialPolygons")) return(TRUE)
    if (is.null(comment(spgeom))) return(TRUE)
    return(comment(spgeom) != "TRUE")
}


createSPComment = function(sppoly,which=NULL,overwrite=TRUE) {
    if (!inherits(sppoly, "SpatialPolygons")) 
        stop("not a SpatialPolygons object")
    if (get_do_poly_check() && notAllComments(sppoly)) { 
      if (is.null(which))
        which = 1:length(sppoly@polygons)
    
      sppoly@polygons[which] = lapply(sppoly@polygons[which], function(p) {
        
        if (!overwrite && !is.null(attr(p, "comment"))) {
            return(p)
        } else if (all(sapply(slot(p, "Polygons"), function(j)
            is.null(slot(j, "coords"))))) {
            comment(p) <- paste(rep(0, length(slot(p, "Polygons"))),
                collapse=" ")
            return(p)
        } else {
            attr(p, "comment") = createPolygonsComment(p)
            return(p)
        }
      })
      comment(sppoly) <- as.character(any(sapply(slot(sppoly, "polygons"),
                function(x) !is.null(comment(x))), na.rm=TRUE))
    }

    return(sppoly)
}

createPolygonsComment = function(poly) {
    if (!is(poly, "Polygons")) 
        stop("not a Polygons object")

    holes = sapply(poly@Polygons, function(x) x@hole)
    if (!any(holes)) {
        comm = rep(0,length(poly@Polygons))
    } else {
        comm = .Call("rgeos_PolyCreateComment", .RGEOS_HANDLE, poly@Polygons, PACKAGE="rgeos")    
    }
    
    return(paste(comm,collapse=" "))
}

getScale <- function() {
    return( mget("scale",.RGEOS_HANDLE)$scale )
}

setScale <- function(scale=100000000) {
    
    maxPreciseValue <- 9007199254740992.0
    
    if(scale > maxPreciseValue){
        stop("Specified scale is larger than maximum allowed")
    }
    
    assign("scale",scale,envir=.RGEOS_HANDLE)
}

checkP4S = function(p4s) {
    
    if ( is.null(p4s) )
        p4s = CRS(as.character(NA))

    if( is.character(p4s))
        p4s = CRS(p4s) 
    
    if (length(p4s) != 1)
        stop("proj4string must be of length 1")
    
    if ( class(p4s) != "CRS") {
        stop("proj4string has invalid class")
    }
    
    return( p4s )
}

translate = function(spgeom) {
    
    rn = row.names(spgeom)
    if (!is.list(rn))
        rn = list(rn)
    
    ids = as.character( unlist( sapply(rn, unique) ) )
    x = .Call("rgeos_double_translate", .RGEOS_HANDLE, spgeom, ids, 0, PACKAGE="rgeos")
    return(x)
}

groupID = function(spgeom, ids) {
    
    if (inherits(spgeom,"SpatialCollections"))
        stop("groupID does not work with SpatialCollections objects")
        
    if (length(row.names(spgeom)) != length(ids)) 
        stop("length of ids does not match number of geometries")
    if (storage.mode(ids) != "character") ids <- as.character(ids)
    
    newids = unique(ids)
    
    if (length(row.names(spgeom)) == 1 || length(row.names(spgeom)) == length(newids) || 
        inherits(spgeom,"SpatialPoints") ) {
        
        row.names(spgeom) <- ids
        return(spgeom)
    }
    
    
    if ( inherits(spgeom,"SpatialLines")  ) {
        
        lineslist = list()
        k=1
        for (curid in newids) {

            if (is.na(curid)) next # RSB 101124
            
            linelist = list()
            l = 1
            for ( i in which(ids == curid) ){
                    
                L = length(spgeom@lines[[i]]@Lines)
                linelist[l:(l+L-1)] = spgeom@lines[[i]]@Lines
                l=l+L
            }
            
            lineslist[[k]] = Lines(linelist, curid)
            k=k+1
        }

        ans = SpatialLines(lineslist,proj4string = spgeom@proj4string)
        
    } else if ( inherits(spgeom,"SpatialPolygons") ) {
        
        polyslist = list()
        k=1
        for (curid in newids) {

            if (is.na(curid)) next # RSB 101124

            comment = c()
            polylist = list()
            l = 1
            for ( i in which(ids == curid) ){
                
                L = length(spgeom@polygons[[i]]@Polygons)
                
                comm = attr(spgeom@polygons[[i]],"comment")
                if (is.null(comm)) comm = rep(0,L)
                else comm = as.integer( strsplit(comm," ")[[1]] )
                comm[comm!=0] = comm[comm!=0] + l-1 
                comment = c(comment, comm)
                    
                polylist[l:(l+L-1)] = spgeom@polygons[[i]]@Polygons
                l=l+L
            }
            
            polyslist[[k]] = Polygons(polylist, curid)
            attr(polyslist[[k]],"comment") = paste(comment, collapse=" ")
            k=k+1
        }

        ans = SpatialPolygons(polyslist,proj4string = spgeom@proj4string)
        
    } else {
        stop("Unknown object class")
    }
    
    return(ans)
}
