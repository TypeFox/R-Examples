# Function to create alpha hull that encompasses x % of occurrences


getDynamicAlphaHull <- function(x, fraction = 0.95, partCount = 3, buff = 10000, coordHeaders = c('Longitude', 'Latitude'), clipToCoast = TRUE, proj = "+proj=longlat +datum=WGS84", verbose = FALSE) {

	if (proj != "+proj=longlat +datum=WGS84") {
		stop("Currently, proj can only be '+proj=longlat +datum=WGS84'.")
	}

	if (ncol(x) == 2) {
		coordHeaders <- c(1,2)
	}

	#reduce to unique coordinates
	x <- x[!duplicated(x[,coordHeaders]), coordHeaders]
	x <- x[complete.cases(x),]
  
	#Alpha hulls cannot be generated if first 3 points are linear. 
	while ((x[1,coordHeaders[1]] == x[2,coordHeaders[1]] & x[2, coordHeaders[1]] == x[3, coordHeaders[1]]) | (x[1,2] == x[2, coordHeaders[2]] & x[2, coordHeaders[2]] == x[3, coordHeaders[2]])) {
		x <- x[sample(1:nrow(x),size = nrow(x)),]
	}

	#create spatialpoints
	x <- SpatialPoints(x, proj4string = CRS(proj))
	x <- remove.duplicates(x)
	if (length(x) < 3) {
		stop('This function requires a minimum of 3 unique coordinates.')
	}

	#create alpha hull and calculate fraction of occurrences that fall within
	#continue until fraction is reached
	alpha = 3
	problem <- FALSE
	if (verbose) {cat('\talpha:', alpha, '\n')}

	hull <- try(alphahull::ahull(data.frame(x),alpha = alpha))
	while ('try-error' %in% class(hull)) {
		if (verbose) {cat('\talpha:', alpha, '\n')}
		alpha <- alpha + 1
		hull <- try(alphahull::ahull(data.frame(x),alpha = alpha))
		if (alpha > 500) {
			problem <- TRUE
			break
		}
	}

	hull <- ah2sp(hull, proj4string = CRS('+proj=longlat +datum=WGS84'))

	if (!is.null(hull)) {
		slot(hull, "polygons") <- lapply(slot(hull, "polygons"), checkPolygonsGEOS2)
	}
 
	while (is.null(hull) | 'try-error' %in% class(hull) | !rgeos::gIsValid(hull, reason = TRUE) %in% c(TRUE, 'Valid Geometry')) {
		alpha <- alpha + 1
		if (verbose) {cat('\talpha:', alpha, '\n')}
		hull <- try(ah2sp(alphahull::ahull(data.frame(x),alpha=alpha), proj4string=CRS('+proj=longlat +datum=WGS84')))
		if (!is.null(hull)) {
			slot(hull, "polygons") <- lapply(slot(hull, "polygons"), checkPolygonsGEOS2)
		}
	}

	#how many points are within hull?
	slot(hull, "polygons") <- lapply(slot(hull, "polygons"), checkPolygonsGEOS2)
	pointWithin <- rgeos::gIntersects(x, hull, byid = TRUE)

	alphaVal <- alpha
	buffered <- FALSE

	while (any(length(hull@polygons[[1]]@Polygons) > partCount, length(which(pointWithin) == TRUE)/length(x) < fraction)) {
	    alpha <- alpha + 1
	    if (verbose) {cat('\talpha:', alpha, '\n')}
	    hull <- try(alphahull::ahull(data.frame(x), alpha = alpha))
	    while ('try-error' %in% class(hull)) {
	      alpha <- alpha + 1
	      hull <- try(alphahull::ahull(data.frame(x),alpha = alpha))
	    }
		hull <- ah2sp(hull, proj4string = CRS('+proj=longlat +datum=WGS84'))
		hull <- sp::spTransform(hull, CRS("+init=epsg:3395"))
		slot(hull, "polygons") <- lapply(slot(hull, "polygons"), checkPolygonsGEOS2)
		hull <- rgeos::gBuffer(hull, width = buff)
		hull <- sp::spTransform(hull, CRS(proj))
		buffered <- TRUE
		pointWithin <- rgeos::gIntersects(x, hull, byid = TRUE)
		alphaVal = alpha
		if (alpha > 500) {
			hull <- rgeos::gConvexHull(x)
			hull <- sp::spTransform(hull, CRS("+init=epsg:3395"))
		    hull <- rgeos::gBuffer(hull, width = buff)
		    hull <- sp::spTransform(hull, CRS(proj))
			buffered <- TRUE
			alphaVal = 'MCH'
			break
		}
	}

	if (!buffered) {
		hull <- sp::spTransform(hull, CRS("+init=epsg:3395"))
		hull <- rgeos::gBuffer(hull,width = buff)
		hull <- sp::spTransform(hull, CRS(proj))
		
	}
  
	if (clipToCoast) {
		# load built-in gshhs dataset
		data(gshhs, envir = environment())
		gshhs <- sp::spTransform(gshhs, CRS(proj4string(hull)))
		hull <- rgeos::gIntersection(hull, gshhs)
	}
  
	return(list(hull, alpha = paste('alpha', alphaVal, sep = '')))
}



#taken from the maptools package
checkPolygonsGEOS2 <- function(obj, properly = TRUE, force = TRUE, useSTRtree = FALSE) {
	if (!is(obj, "Polygons")) 
	    stop("not a Polygons object")
	comm <- try(rgeos::createPolygonsComment(obj), silent = TRUE)
	if (class(comm) != "try-error" && !force) {
	    comment(obj) <- comm
	    return(obj)
	}
	pls <- slot(obj, "Polygons")
	IDs <- slot(obj, "ID")
	n <- length(pls)
	if (n < 1) 
	    stop("Polygon list of zero length")
	uniqs <- rep(TRUE, n)
	if (n > 1) {
	    if (useSTRtree) 
	        tree1 <- rgeos::gUnarySTRtreeQuery(obj)
	    SP <- SpatialPolygons(lapply(1:n, function(i) Polygons(list(pls[[i]]), ID = i)))
	    for (i in 1:(n - 1)) {
	        if (useSTRtree) {
	            if (!is.null(tree1[[i]])) {
	              res <- try(rgeos::gEquals(SP[i, ], SP[tree1[[i]],], byid = TRUE), silent = TRUE)
	              if (class(res) == "try-error") {
	                warning("Polygons object ", IDs, ", Polygon ", i, ": ", res)
	                next
	              }
	              if (any(res)) {
	                uniqs[as.integer(rownames(res)[res])] <- FALSE
	              }
				}
	        }
	        else {
	            res <- try(rgeos::gEquals(SP[i, ], SP[uniqs,], byid = TRUE), silent = TRUE)
	            if (class(res) == "try-error") {
	              warning("Polygons object ", IDs, ", Polygon ", i, ": ", res)
	              next
	            }
	            res[i] <- FALSE
	            if (any(res)) {
	              wres <- which(res)
	              uniqs[wres[wres > i]] <- FALSE
	            }
	        }
	    }
	}
	if (any(!uniqs)) 
	    warning(paste("Duplicate Polygon objects dropped:", paste(wres, collapse = " ")))
	pls <- pls[uniqs]
	n <- length(pls)
	if (n < 1) 
	    stop("Polygon list of zero length")
	if (n == 1) {
	    oobj <- Polygons(pls, ID = IDs)
	    comment(oobj) <- rgeos::createPolygonsComment(oobj)
	    return(oobj)
	}
	areas <- sapply(pls, slot, "area")
	pls <- pls[order(areas, decreasing = TRUE)]
	oholes <- sapply(pls, function(x) slot(x, "hole"))
	holes <- rep(FALSE, n)
	SP <- SpatialPolygons(lapply(1:n, function(i) Polygons(list(pls[[i]]), ID = i)))
	if (useSTRtree) 
	    tree2 <- rgeos::gUnarySTRtreeQuery(SP)
	for (i in 1:(n - 1)) {
	    if (useSTRtree) {
	        if (!is.null(tree2[[i]])) {
	            if (properly) 
	              res <- rgeos::gContainsProperly(SP[i, ], SP[tree2[[i]], ], byid = TRUE)
	            else res <- rgeos::gContains(SP[i, ], SP[tree2[[i]], ], byid = TRUE)
	        }
	        else {
	            res <- FALSE
	        }
	    }
	    else {
	        if (properly) 
	            res <- rgeos::gContainsProperly(SP[i, ], SP[-(1:i), ], byid = TRUE)
	        else res <- rgeos::gContains(SP[i, ], SP[-(1:i), ], byid = TRUE)
	    }
	    wres <- which(res)
	    if (length(wres) > 0L) {
	        nres <- as.integer(rownames(res))
	        holes[nres[wres]] <- !holes[nres[wres]]
	    }
	}
	for (i in 1:n) {
	    if (oholes[i] != holes[i]) 
	        pls[[i]] <- Polygon(slot(pls[[i]], "coords"), hole = holes[i])
	}
	oobj <- Polygons(pls, ID = IDs)
	comment(oobj) <- rgeos::createPolygonsComment(oobj)
	oobj
}






