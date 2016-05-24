# Author: Robert J. Hijmans
# Date : December 2009
# Version 0.1
# Licence GPL v3


.randomCellsLonLat <- function(r, n) {
# sampling cells with weights based on their acutal area, to avoid sampling 
# too many cells when going towards the poles
# as suggested by Jane Elith
	y <- yFromRow(r, 1:nrow(r)) # get y coordinates for all rows
	dx <- pointDistance(cbind(0, y), cbind(xres(r), y), longlat=TRUE)  # size of each cell in longitude directoin
	dx <- dx / max(dx) # standardize to 1

	row <- sample.int(nrow(r), n, replace = TRUE, prob = dx) # sample from rows using weights
	col <- sample.int(ncol(r), n, replace = TRUE) # sample from cols, not using weights
	rowcol <- unique(cbind(row, col))

	maxrow <- pmax(1, round(dx * ncol(r)))
	cellsel <- matrix(nrow=0, ncol=2)
	for (i in unique(rowcol[,1])) {
		a <- subset(rowcol, rowcol[,1]==i)
		if (nrow(a) > maxrow[i]) { a <- a[1:maxrow[i]] }
		cellsel <- rbind(cellsel, a)
	}

	cells <- cellFromRowCol(r, cellsel[,1], cellsel[,2])
	return(cells)
}



..randomCellsLonLat2 <- function(r, n) {
# only cells that are not NA
	cells <- which(! is.na(getValues(r)) )
# which rows are that?
	rows <- rowFromCell(r, cells)
# what is the latitude?
	y <- yFromRow(r, rows)
# what is the 'width' of a cell?
	dx <- pointDistance(cbind(0, y), cbind(xres(r), y), longlat=TRUE)  

	cells <- sample(cells, n, prob=dx)
	return(cells)
}



randomPoints <- function(mask, n, p, ext=NULL, extf=1.1, excludep=TRUE, prob=FALSE, cellnumbers=FALSE, tryf=3, warn=2, lonlatCorrection=TRUE) {
	
	if (nlayers(mask) > 1) { 
		mask <- raster(mask, 1)	
	}
	
	tryf <- max(round(tryf[1]), 1)
	
	if (missing(p)) { 
		excludep <- FALSE
	} else {
		if (inherits(p, 'SpatialPoints')) {
			p <- coordinates(p)
		}
	}
	
	if (class(ext)=='character') {
		if (! ext %in% c('points')) { 
			stop("if ext is a character variable it should be 'points'") 
		} else if (missing(p)) { 
			warning("if p is missing, 'ext=points' is meaningless") 
			ext <- extent(mask)  
		} else {
			ext <- extent(min(p[,1]), max(p[,1]), min(p[,2]), max(p[,2]))
		}
	} 

	if (! is.null(ext)) {
		ext <- extent(ext)
		ext <- ext * extf
		ext <- intersect(ext, extent(mask))
		mask2 <- crop(raster(mask), ext)
	}  else {
		mask2 <- raster(mask)
	}
	
	if (n > ncell(mask2)) {
		n <- ncell(mask2)
		if (warn>0) { warning('changed n to ncell of the mask (extent)') }
	}
	
	nn <- n * tryf
	nn <- max(nn, 10)

	if (prob) {
		stopifnot(hasValues(mask))
		cells <- crop(mask, mask2)
		cells <- try( stats::na.omit(cbind(1:ncell(cells), getValues(cells))))
		if (class(cells) == 'try-error') {
			stop("the raster is too large to be used with 'prob=TRUE'")
		}
		prob <- cells[,2]
		cells <- cells[,1]
		if (couldBeLonLat(mask)) {
			rows <- rowFromCell(mask2, cells)
			y <- yFromRow(mask2, rows)
			dx <- pointDistance(cbind(0, y), cbind(xres(mask2), y), longlat=TRUE)  
			prob <- prob * dx
		}

		cells <- sample(cells, nn, prob=prob)
		xy <- xyFromCell(mask2, cells)
		cells <- cellFromXY(mask, xy)
		
	} else 	if (canProcessInMemory(mask2)) {
	
		cells <- crop(mask, mask2)
		if (hasValues(cells)) {
			cells <- which(! is.na(getValues(cells)) )
		} else {
			cells <- 1:ncell(cells)
		}
		nn <- min(length(cells), nn)

		if (lonlatCorrection & couldBeLonLat(mask)) {
			# which rows are that?
			rows <- rowFromCell(mask2, cells)
			# what is the latitude?
			y <- yFromRow(mask2, rows)
			# what is the 'width' of a cell?
			dx <- pointDistance(cbind(0, y), cbind(xres(mask2), y), longlat=TRUE)  
			cells <- sample(cells, nn, prob=dx)
	
		} else {
			cells <- sample(cells, nn)
		}
		xy <- xyFromCell(mask2, cells)
		cells <- cellFromXY(mask, xy)
	
	} else {
	
		nn <- min(ncell(mask2), nn)
		if (couldBeLonLat(mask2)) {
	
			cells <- .randomCellsLonLat(mask2, nn)
			
		} else {
			if (nn >= ncell(mask2)) {
				cells <- 1:ncell(mask2)
			} else {
				cells <- sampleInt(ncell(mask2), nn)
			}
		
		}
		xy <- xyFromCell(mask2, cells)
		cells <- cellFromXY(mask, xy)
		
		if (hasValues(mask)) {
			
			vals <- cbind(cells, extract(mask, cells))
			cells <- stats::na.omit(vals)[,1]
		}
	}
		
	if (excludep) {	
		pcells <- cellFromXY(mask, p)
		cells <- cells[ ! (cells %in% pcells) ] 	
	}

	if (length(cells) > n) { 
		
		cells <- sample(cells, n) 
		
	} else if (length(cells) < n) {
	
		frac <- length(cells) / n
		if (frac < 0.1) {
			stop("generated random points = ", frac," times requested number; Use a higher value for tryf" )
		}
		if (frac < 0.5  & warn==1) {
			warning("generated random points = ", frac," times requested number; Use a higher value for tryf" )
		} else if (warn > 1) {
			warning("generated random points = ", frac," times requested number")
		}
	}
	
	if (cellnumbers) {
		return(cells)
	} else {
		return(xyFromCell(mask, cells))
	}
}



#	if (canProcessInMemory(mask, 2)) {
#		if (dataContent(mask) != 'all') { mask <- readAll(mask) }
#		if (e < extent(mask)) { mask <- crop(mask, e) }
#		if (dataContent(mask) == 'all') {
#			cells <- stats::na.omit(cbind(1:ncell(r), values(mask)))[,1]
#		} else {
#			cells <- 1:ncell(mask)
#		}
#		if (excludep) {	
#			pcells <- cellFromXY(mask, p)
#			cells <- cells[!(cells%in%pcells)] 	}
#		}
#		if (n > length(cells)) {
#			warning('there are only ',length(cells),' cells available')	
#		} else {
#			cells <- sample(cells, n)
#		}

