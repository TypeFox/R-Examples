
# taken from the raster package
#   the function isn't exported from raster so it must be replicated here
# author Robert Hijmans
# June 2010
# version 1.0
# license GPL3


.compareCRS <- function(x, y, unknown=FALSE, verbatim=FALSE, verbose=FALSE) {
	
	x <- tolower(projection(x))
	y <- tolower(projection(y))
	
	step1 <- function(z, verbatim) {
		z <- gsub(' ', '', z)
		if (!verbatim) {
			z <- unlist( strsplit(z, '+', fixed=TRUE) )[-1]
			z <- do.call(rbind, strsplit(z, '='))
		}
		z
	}
	
	if (verbatim) {
		return(x==y)
	}
	
	x <- step1(x, verbatim)
	y <- step1(y, verbatim)
	
	if (length(x) == 0 | length(y) == 0) {
		if (unknown) {
			return(TRUE)
		} else {
			if (verbose) {
				cat('Unknown CRS\n')
			}
			return(FALSE) 
		}
	}
	x <- x[x[,1] != 'towgs84', , drop=FALSE]
	x <- x[which(x[,1] %in% y[,1]), ,drop=FALSE]
	y <- y[which(y[,1] %in% x[,1]), ,drop=FALSE]
	x <- x[order(x[,1]), ,drop=FALSE]
	y <- y[order(y[,1]), ,drop=FALSE]
	i <- x[,2] == y[,2]
	
	if (! all(i)) {
		if (verbose) {
			i <- which(!i)
			for (j in i) {
				cat('+',x[j,1], ':  ', x[j,2],' != ', y[j,2], '\n', sep='') 
			}
		}
		return(FALSE)
	}
	return(TRUE)
}

cellsBuffer = function(cells, buffer) {
	cells = squareRaster(cells)
  
	buffer =  ceiling(buffer/xres(cells))
	
	cellsInla = raster::extend(cells, c(buffer, buffer))
	values(cellsInla ) =  
			c(t(matrix(seq(1,ncell(cellsInla)), 
									nrow=nrow(cellsInla), ncol=ncol(cellsInla))))
	names(cellsInla) = "space"
	
	
	result = list(small=crop(cellsInla, cells), big=cellsInla)
}

