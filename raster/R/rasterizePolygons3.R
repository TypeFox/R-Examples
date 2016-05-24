# Author: Robert J. Hijmans
# Date : July 2011
# Version 1.0
# Licence GPL v3


.checkFields <- function(putvals) {
	for (i in 1:ncol(putvals)) {
		if (class(putvals[,i]) == 'factor') {
			ptv <- try( as.numeric(as.character(putvals[,i])) )
			if (class(ptv) == 'try-error') {
				putvals[,i] <- as.numeric(putvals[,i])
			} else {
				putvals[,i] <- ptv
			}
		} else if (class(putvals) == 'character') {
			ptv <- try( as.numeric(as.character(putvals[,i])) )
			if (class(ptv) == 'try-error') {
				stop('field: ', colnames(putvals)[i], ' cannot be converted to a number')
			} else {
				putvals <- ptv
			}
		}
	}
	putvals
}


.p3r <- function(p, x, field=0, background=NA, mask=FALSE, update=FALSE, filename="", ...) {

	filename <- trim(filename)
	if (mask | update) {
		if (mask & update) {
			stop('you cannot use mask=TRUE & update=TRUE at the same time')
		}
		oldx <- x
		
	}
	x <- raster(x)
	projp <- projection(p)
	if (! compareCRS(projp, x)) {
	#	warning('crs or raster and polygons do not match')
	}
	if (! is.na(projp)) {
		projection(x) = projp
	} 
	
	spbb <- bbox(p)
	rsbb <- bbox(x)
	if (spbb[1,1] >= rsbb[1,2] | spbb[1,2] <= rsbb[1,1] | spbb[2,1] >= rsbb[2,2] | spbb[2,2] <= rsbb[2,1]) {
		stop('polygon and raster have no overlapping areas')
	}
	npol <- length(p@polygons)
	
	if (! is.numeric(field) ) {
		putvals <- .checkFields( p@data[, field, drop=FALSE] )
	} else if (length(field) > 1) { 
		if (length(field) == npol) {
			putvals <- field
		} else {
			if (inherits(p, 'SpatialPolygonsDataFrame')) {
				putvals <- .checkFields( p@data[, field, drop=FALSE] )
			} else {
				stop('incorrect value for "field"')
			}
		}
	} else if (mask) {
		putvals <- rep(1, length=npol)	
	} else if (field < 0) {
		putvals <- rep(1, length=npol)	
	} else if (field == 0) {
		putvals <- as.integer(1:npol)
	} else {
		putvals <- .checkFields(p@data[, field, drop=FALSE])
	}
	if (is.vector(putvals)) {
		putvals <- matrix(putvals, ncol=1)
	} else {
		putvals <- as.matrix(putvals)
	}

	nl <- NCOL(putvals)
	if (nl > 1) {
		x <- brick(x, nl=nl)
	}
	
	npol <- length(p@polygons)
	polinfo <- matrix(NA, nrow=npol * 3, ncol=7)
	addpol <- matrix(NA, nrow=500, ncol=7)
	polx <- poly <- vector(length=npol * 3, mode='list')
	cnt <- 0
	
	for (i in 1:npol) {
		nsubpol <- length(p@polygons[[i]]@Polygons)
		for (j in 1:nsubpol) {
			cnt <- cnt + 1
			if (cnt > dim(polinfo)[1]) { 
				polinfo <- rbind(polinfo, addpol)  
			}
			polinfo[cnt, 1] <- cnt
			rg <- range(p@polygons[[i]]@Polygons[[j]]@coords[,1])
			polinfo[cnt, 2] <- rg[1]
			polinfo[cnt, 3] <- rg[2]
			rg <- range(p@polygons[[i]]@Polygons[[j]]@coords[,2])
			polinfo[cnt, 4] <- rg[1]
			polinfo[cnt, 5] <- rg[2]
			polinfo[cnt, 6] <- p@polygons[[i]]@Polygons[[j]]@hole 
			polinfo[cnt, 7] <- i
			
			polx[[cnt]] <- p@polygons[[i]]@Polygons[[j]]@coords[,1]
			poly[[cnt]] <- p@polygons[[i]]@Polygons[[j]]@coords[,2]
		}
	}
	#rm(p)
	polinfo <- subset(polinfo, polinfo[,1] <= cnt, drop=FALSE)
	
	polx <- polx[1:cnt]
	poly <- poly[1:cnt]
	
	message('Found', npol, 'region(s) and', cnt, 'polygon(s)') 
		
	if (!canProcessInMemory(x)) {
		if (filename == "") {
			filename <- rasterTmpFile()
		}
	}
	if (filename == "") {
		v <- matrix(NA, ncol=nlayers(x), nrow=ncell(x))
	} else {
		x <- writeStart(x, filename=filename, ...)
	}

	tr <- blockSize(x, n=2*nlayers(x))
	pb <- pbCreate(tr$n, label='rasterize', ...)

	rx <- c(xmin(x), xmax(x))
	
	for (i in 1:tr$n) {
		
		cells <- cellFromRowCol(x, tr$row[i], 1) : cellFromRowCol(x, tr$row[i]+(tr$nrows[i]-1), ncol(x))
		xy <- xyFromCell(x, cells)
		ry <- range(xy[,2])
		subpol <- subset(polinfo, !(polinfo[,2] > rx[2] | polinfo[,3] < rx[1] | polinfo[,4] > ry[2] | polinfo[,5] < ry[1] ), drop=FALSE)

		if (nrow(subpol) > 0) { 		
			rrv <- rep(NA, nrow(xy))
				
			px <- polx[subpol[,1]]
			py <- poly[subpol[,1]]
			p <- xy[,1] >= min(subpol[,2]) &  xy[,1] <= max(subpol[,3]) &  xy[,2] >= min(subpol[,4]) & xy[,2] <= max(subpol[,5])

			rrv[p] <- .Call('point_in_polygon2', xy[p,1], xy[p,2], px, py, as.integer(subpol[,7]), as.integer(subpol[,6]), PACKAGE='raster')
		
			if (!is.na(background)) { 
				rrv[is.na(rrv)] <- background
			}
		
			if (mask) {
				vals <- getValues(oldx, tr$row[i], tr$nrows[i])
				if (nl == 1) {
					vals <- matrix(vals, ncol=1)
				}
				vals[is.na(rrv), ] <- NA
			} else if (update) {
				vals <- getValues(oldx, tr$row[i], tr$nrows[i])
				if (nl == 1) {
					vals <- matrix(vals, ncol=1)
				}
				vals[is.na(rrv), ] <- putvals[!is.na(rrv), ]
			} else {
				vals <- putvals[rrv, ]
			}
			
			if (filename == "") {
				v[cells,] <- vals
			} else {
				x <- writeValues(x, vals, tr$row[i])
			}
		} else {
			if (filename != "") {
				vals <- matrix(NA, nrow=length(cells), ncol=nl)
				x <- writeValues(x, vals, tr$row[i])
			}
		}
		pbStep(pb, i)
	}
	pbClose(pb)

	if (filename == "") {
		x <- setValues(x, v)
	} else {
		x <- writeStop(x)
	}
	return(x)
}

#e = .p3r(p, r)

 
 