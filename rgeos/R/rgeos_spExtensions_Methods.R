#########################
###
### SpatialRings
###
#########################

Ring <- function(coords,ID=as.character(NA)) {
    if (ncol(coords) != 2) 
        stop("coords must be a two-column matrix")
    
    n = nrow(coords)
    area2 = sum( (coords[-1,1]-coords[-n,1])*(coords[-1,2]+coords[-n,2]) )
    if (area2 < 0) {
        # if area2 is negative coordinates are ccw, reverse them
        coords[,1] = rev(coords[,1])
        coords[,2] = rev(coords[,2])
    }
    
    coords <- coordinates(coords)
    new("Ring", coords = coords, ID = ID)
}

.bboxSRs <- function(lst) {
    bb = sapply(lst, bbox)
    res = matrix(c(min(bb[1,]), min(bb[2,]), max(bb[3,]), max(bb[4,])), 2, 2)
    dimnames(res) = list(c("x", "y"), c("min", "max"))
    return(res)
}

SpatialRings <- function(RingList, proj4string=CRS(as.character(NA))) {
    if (any(sapply(RingList, function(x) !is(x, "Ring")))) 
        stop("Ring list not exclusively filled with Ring objects")
    Sp <- new("Spatial", bbox = .bboxSRs(RingList), proj4string=proj4string)
    res <- new("SpatialRings", Sp, rings=RingList)
    return(res)
}


bbox.Ring <- function(obj) {
    rx <- range(obj@coords[,1])
    ry <- range(obj@coords[,2])
    res = rbind(r1 = rx, r2 = ry)
    dimnames(res)[[2]] <- c("min", "max")
    return(res)
}
setMethod("bbox", "Ring", bbox.Ring)


plotSpatialRings <- function(SR, xlim = NULL, ylim = NULL,
                             col = 1, lwd = 1, lty=1, add = FALSE, axes = FALSE, ..., 
                             setParUsrBB=FALSE) {

    if (!add) 
        plot(as(SR, "Spatial"), xlim = xlim, ylim = ylim, axes = axes, ..., setParUsrBB=setParUsrBB)
    

    lst <- SR@rings
    if (length(col) != length(lst)) 
        col <- rep(col[1], length(lst))
    if (length(lwd) != length(lst)) 
        lwd <- rep(lwd[1], length(lst))
    if (length(lty) != length(lst)) 
        lty <- rep(lty[1], length(lst))

    for (i in seq(along=lst)) {
        crds <- coordinates(lst[[i]])
        lines(crds, col = col[i], lwd = lwd[i], lty = lty[i], ...)
    }
}

setMethod("plot", 
          signature(x = "SpatialRings", y = "missing"), 
          function(x, y, ...) plotSpatialRings(x, ...) )

setMethod("coordinates", "Ring", function(obj) obj@coords)
setMethod("coordinates", "SpatialRings", function(obj) lapply(obj@rings, coordinates))


#if (!isGeneric("lines"))
#    setGeneric("lines", function(x, y, ...)
#        standardGeneric("lines"))

#setMethod("lines", "Ring", function(x, y = NULL, ...) invisible(lines(coordinates(x), ...)))
#setMethod("lines", "SpatialRings", function(x, y = NULL, ...) {
#    f = function(x, ...) lines(x, ...)
#    invisible(lapply(x@rings, f, ...))
#})

if (!isGeneric("row.names"))
	setGeneric("row.names", function(x) standardGeneric("row.names"))
if (!isGeneric("row.names<-")) 
	setGeneric("row.names<-", function(x, value) standardGeneric("row.names<-"))

setMethod("row.names", "SpatialRings", function(x) sapply(slot(x, "rings"), slot, "ID"))
setReplaceMethod("row.names", signature(x = "SpatialRings", value = "character"),
                 function(x, value) spChFIDs(x, value))

setMethod("[", "SpatialRings", 
    function(x, i, j, ..., drop = TRUE) {
        if (any(is.na(i))) stop("NAs not permitted in row index")
        
        if (is.logical(i)) {
            if (length(i) == 1 && i) {
                i = 1:length(x@rings)
            } else {
                i <- which(i)
            }
        } else if (is.character(i)) {
            i <- match(i, row.names(x))
        }
        
        x@rings = x@rings[i]
        x@bbox = .bboxSRs(x@rings)
        return(x)
    }
)


setMethod("coordnames", signature(x = "SpatialRings"), function(x) coordnames(x@rings[[1]]))
setMethod("coordnames", signature(x = "Ring"), function(x) dimnames(coordinates(x))[[2]])

setReplaceMethod("coordnames", 
                 signature(x = "SpatialRings", value = "character"),
                 function(x, value) {
                     dimnames(x@bbox)[[1]] = value
                     for (i in seq(along = x@rings))
                        coordnames(x@rings[[i]]) = value
                     return(x)
                 }
)
setReplaceMethod("coordnames",signature(x = "Ring", value = "character"),
                 function(x, value) {
                     dimnames(x@coords)[[2]] = value
                     return(x)
                 }
)



chFIDsSpatialRings <- function(obj, x) {
    nl <- length(slot(obj, "rings"))
    if (length(x) != nl) stop("lengths differ")
    if (length(x) > length(unique(x))) stop("duplicate IDs")

    rings <- slot(obj, "rings")
    for (i in 1:nl) obj@rings[[i]]@ID = x[i]
    
    return(obj)
}
setMethod("spChFIDs", signature(obj="SpatialRings", x="character"), chFIDsSpatialRings)

setAs("Ring", 
      "SpatialPoints",
      function(from) { 
          SpatialPoints(do.call("rbind", coordinates(from)))
      }
)

setAs("SpatialRings", 
      "SpatialPoints", 
      function(from) { 
          SpatialPoints(do.call("rbind", 
                        lapply(from@rings, function(x) as(x, "SpatialPoints"))),
                        CRS(proj4string(from)))
      }
)


#########################
###
### SpatialRingsDataFrame
###
#########################

SpatialRingsDataFrame = function(sr, data, match.ID = TRUE) {
    if (match.ID) {
        sr_IDs <- row.names(sr)
        data_IDs <- row.names(data)
        mtch <- match(sr_IDs, data_IDs)
        if (any(is.na(mtch)))
            stop("row.names of data and Rings IDs do not match")
        if (length(unique(mtch)) != length(sr_IDs))
            stop("row.names of data and Rings IDs do not match")
        data <- data[mtch, , drop=FALSE]
    }
    if (nrow(data) != length(sr@rings))
        stop("length of data.frame does not match number of Ring elements")
    
    return( new("SpatialRingsDataFrame", sr, data = data) )
}

setAs("SpatialRingsDataFrame","SpatialRings", function(from) SpatialRings(from@rings))
setAs("SpatialRingsDataFrame", "data.frame", function(from) from@data)

setMethod("names","SpatialRingsDataFrame", function(x) names(x@data))    
setReplaceMethod("names", signature(x = "SpatialRingsDataFrame", value = "character"),
                 function(x, value) { names(x@data)<-value; x })

setMethod("row.names","SpatialRingsDataFrame", function(x) sapply(slot(x, "rings"), slot, "ID"))    
setReplaceMethod("row.names", signature(x = "SpatialRingsDataFrame", value = "character"),
              function(x, value) spChFIDs(x, value))


setMethod("[", c("SpatialRingsDataFrame", "ANY", "ANY"), 
    function(x, i, j, ... , drop = TRUE) {
        missing.i = missing(i)
        missing.j = missing(j)
        nargs = nargs() # e.g., a[3,] gives 2 for nargs, a[3] gives 1.
        
        if (missing.i && missing.j) {
            i = TRUE
            j = TRUE
        } else if (missing.j && !missing.i) {
            if (nargs == 2) {
                j = i
                i = TRUE
            } else {
                j = TRUE
            }
        } else if (missing.i && !missing.j) {
            i = TRUE
        }
        
        if (is.matrix(i)) stop("matrix argument not supported in SpatialRingsDataFrame selection")
        if (is.logical(i)) {
            if (length(i) == 1 && i)
                i = 1:length(x@rings)
            else
                i = which(i)
        } else if (is.character(i)) {
                i = match(i, row.names(x))
        }
        
        if (any(is.na(i))) stop("NAs not permitted in row index")
        
        x@rings = x@rings[i]
        x@data = x@data[i, j, ..., drop = FALSE]
        x@bbox = .bboxSRs(x@rings)
        
        return(x)
    }
)

#setMethod("lines","SpatialRingsDataFrame", function(x, y = NULL, ...) lines(as(x, "SpatialRings"), ...))    
setMethod("dim","SpatialRingsDataFrame", function(x) dim(x@data))


chFIDsSpatialRingsDataFrame <- function(obj, x) {
    SR <- as(obj, "SpatialRings")
    SRx <- spChFIDs(SR, x)
    df <- as(obj, "data.frame")
    row.names(df) <- sapply(slot(SRx, "rings"), function(x) slot(x, "ID"))
    SpatialRingsDataFrame(SRx, data=df)
}

setMethod("spChFIDs", signature(obj="SpatialRingsDataFrame", x="character"), chFIDsSpatialRingsDataFrame)



#########################
###
### SpatialCollections
###
#########################


SpatialCollections <- function( points = NULL, lines = NULL,
								rings = NULL, polygons = NULL,
								plotOrder = c(4,3,2,1),
								proj4string=CRS(as.character(NA))) {

	plotOrder=as.integer(plotOrder)
	stopifnot(is.integer(plotOrder))
	stopifnot(length(plotOrder) == 4)
	
	stopifnot(inherits(points,"SpatialPoints") | is.null(points))
	stopifnot(inherits(lines,"SpatialLines") | is.null(lines))
	stopifnot(inherits(rings,"SpatialRings") | is.null(rings))
	stopifnot(inherits(polygons,"SpatialPolygons") | is.null(polygons))
	stopifnot(is(proj4string, "CRS"))
	
	bb = c()
	if (!is.null(points)) bb = rbind(bb,c(bbox(points)))
	if (!is.null(lines)) bb = rbind(bb,c(bbox(lines)))
	if (!is.null(rings)) bb = rbind(bb,c(bbox(rings)))
	if (!is.null(polygons)) bb = rbind(bb,c(bbox(polygons)))
	
	if (length(bb) == 0) {
		bbox = matrix( rep(NA, 2), 2, 2)
	} else {
		bbox = matrix(c(min(bb[,1]), min(bb[,2]), max(bb[,3]), max(bb[,4])), 2, 2)
	}
	dimnames(bbox) = list(c("x", "y"), c("min", "max"))
	
	
	
	Sp <- new("Spatial", bbox = bbox, proj4string=proj4string)	
    res <- new("SpatialCollections", Sp, plotOrder=plotOrder, pointobj=points, 
				lineobj=lines, ringobj=rings, polyobj=polygons)
	#validObject(res)

    return(res)
}


plotSpatialCollections <- function(SC, 
								   pointopt = list(), lineopt = list(),
								   ringopt = list(), polyopt = list(),
								   pch = 3, cex = 1, bg = 1,
								   col = 1, lwd = 1, lty=1,
								   border = par("fg"), xpd = NULL, 
								   density = NULL, angle = 45, pbg=NULL,
								   xlim = NULL, ylim = NULL,
								   add = FALSE, axes = FALSE, ...,
	  							   setParUsrBB=FALSE) {
		
	if (!add) {
		plot(as(SC,"Spatial"), xlim=xlim, ylim=ylim, axes=axes, 
			 ..., setParUsrBB=setParUsrBB)
	}
	add=TRUE
	
	call = match.call()
	if (is.null(call$col) & is.null(polyopt[["col"]]))
		polyopt[["col"]] = NA
		
	for (i in order(SC@plotOrder)) {
		
		if (i == 1 & !is.null(SC@pointobj)) { # plot points
			ptcol = col
			ptlwd = lwd
			
			if (!is.null(pointopt[["pch"]])) ptpch = pointopt[["pch"]]
			if (!is.null(pointopt[["cex"]])) ptcex = pointopt[["cex"]]
			if (!is.null(pointopt[["col"]])) ptcol = pointopt[["col"]]
			if (!is.null(pointopt[["lwd"]])) ptlwd = pointopt[["lwd"]]
			if (!is.null(pointopt[["bg"]]))  ptbg  = pointopt[["bg"]]
			
			plot(SC@pointobj, pch=pch, axes=axes, add=add, xlim=xlim, ylim=ylim, 
				 cex=cex, col=ptcol, lwd=ptlwd, bg=bg, ..., setParUsrBB=setParUsrBB)
			
		} else if (i == 2 & !is.null(SC@lineobj)) { #plot lines
			lcol = col
			llwd = lwd
			llty = lty
			
			if (!is.null(lineopt[["col"]])) lcol = lineopt[["col"]]
			if (!is.null(lineopt[["lwd"]])) llwd = lineopt[["lwd"]]
			if (!is.null(lineopt[["lty"]])) llty = lineopt[["lty"]]
			
			plot(SC@lineobj, xlim=xlim, ylim=ylim, col=lcol, lwd=llwd, lty=llty,
				 add=add, axes=axes, ..., setParUsrBB=setParUsrBB)
			
		} else if (i == 3 & !is.null(SC@ringobj)) { #plot rings
			rcol = col
			rlwd = lwd
			rlty = lty
			
			if (!is.null(ringopt[["col"]])) rcol = ringopt[["col"]]
			if (!is.null(ringopt[["lwd"]])) rlwd = ringopt[["lwd"]]
			if (!is.null(ringopt[["lty"]])) rlty = ringopt[["lty"]]
			
			plot(SC@ringobj, xlim=xlim, ylim=ylim, col=rcol, lwd=rlwd, lty=rlty,
				 add=add, axes=axes, ..., setParUsrBB=setParUsrBB)
				
		} else if (i == 4 & !is.null(SC@polyobj)) { #plot polygons
			pcol = col
			if (!is.null(polyopt[["col"]])) pcol = polyopt[["col"]]
			
			plot(SC@polyobj, col=pcol, border=border, add=add, xlim=xlim, ylim=ylim,
				 xpd=xpd, density=density, angle=angle, pbg=pbg, axes=axes, ...,
				 setParUsrBB=setParUsrBB)
		}
	}
}
		
setMethod("plot", signature(x = "SpatialCollections", y = "missing"),
	function(x, y, ...) plotSpatialCollections(x, ...))
	
setMethod("row.names", "SpatialCollections", function(x) {
	
	ans = list()
	if (!is.null(x@pointobj)) ans$points = row.names(x@pointobj)
	if (!is.null(x@lineobj)) ans$lines = row.names(x@lineobj)
	if (!is.null(x@ringobj)) ans$rings = row.names(x@ringobj)
	if (!is.null(x@polyobj)) ans$polygons = row.names(x@polyobj)
	
	return(ans)
})
