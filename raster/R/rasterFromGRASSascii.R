# Author: Robert J. Hijmans
# Date : October 2009
# Version 0.9
# Licence GPL v3

# Changes for GRASS by Cornel M. Pop <cornel@popgeology.com> 


.rasterFromGRASSASCIIFile <- function(filename, crs=NA, ...) {

  offset <- 9 # Maximum real offset (GRASS = 6 required + 3 opts; ESRI = 5 required + 1 opt)
  #offset <- as.integer(offset)
  
	splitasc <- function(s) {
		s <- trim(s)
		spl <- unlist(strsplit(s, ''), use.names = FALSE)
		pos <- which(spl==' ')[1]
		first <- substr(s, 1, (pos-1))
		second <- substr(s, (pos+1), .nchar(s))
		return(trim(c(first, second)))
	}
	
	filename <- trim(filename)
    if (!file.exists(filename)) { stop(paste(filename, " does not exist")) }
	con <- file(filename, "rt")
	lines <- readLines(con, n=offset)
	close(con)
	ini <- lapply(lines, splitasc)
  ini <- matrix(unlist(ini, use.names = FALSE), ncol=2, byrow=TRUE)

  # Data may begin with a noval, which may not be numeric (i.e. NaN, NA)
  noval_kw <- c("NULL", "NODATA", "NODATA_VALUE") # Note that NODATA is not in the ESRI specs - using here for compatibility
  noval <- ini[which(toupper(ini[,1]) %in% noval_kw == T), 2]
  if(length(noval) == 1){
    offset <- length(which(is.na(suppressWarnings(try(as.numeric(ini[which(ini[,1] != noval),1])))))) # True offset. Rest is overshot (data)
  } else {
    offset <- length(which(is.na(suppressWarnings(try(as.numeric(ini[,1]))))))
  }
  stopifnot(offset > 2) # Note: This is not necessary.
  ini <- ini[1:offset,] # Keep only real header elements
	ini[,1] = toupper(ini[,1]) 
  
	nodataval <- xn <- yn <- d <- nr <- nc <- xc <- yc <- NA
  grass_mandatory <- c("NORTH:", "SOUTH:", "EAST:", "WEST:", "ROWS:", "COLS:")
  if (!FALSE %in% (grass_mandatory %in% ini[,1])){
    # GRASS ASCII - note: this is a bit different from how ESRI files are
    # handled, in that we check for all mandatory header keys, and don't process
    # the file if they are missing because, well, they should not be missing!
    nr <- as.integer(ini[which(ini[,1] == "ROWS:"),2])
    nc <- as.integer(ini[which(ini[,1] == "COLS:"),2])
    xn <- as.numeric(ini[which(ini[,1] == "WEST:"),2])
    yn <- as.numeric(ini[which(ini[,1] == "SOUTH:"),2])
    xx <- as.numeric(ini[which(ini[,1] == "EAST:"),2])
    yx <- as.numeric(ini[which(ini[,1] == "NORTH:"),2])
    
    # Optional fields:
    if ("NULL:" %in% ini[,1]){
      nodataval <- as.numeric(ini[which(ini[,1] == "NULL:"),2])
    } else {
      warning('"NULL" not detected. Setting it to -Inf\n  You can set it to another value with function "NAvalue"')
      nodataval <- -Inf
    }
    if ("TYPE:" %in% ini[,1]){
      warning('Optional "TYPE" parameter detected, but ignoring and using the default data type.')
    }
    if ("MULTIPLIER:" %in% ini[,1]){multiplier <- as.numeric(ini[which(ini[,1] == "MULTIPLIER:"),2])}

  } else if ("NCOLS" %in% ini[,1]){
	  for (i in 1:nrow(ini)) {
		  if (ini[i,1] == "NCOLS") { nc <- as.integer(ini[i,2])
		  } else if (ini[i,1] == "NROWS") { nr <- as.integer(ini[i,2])
		  } else if (ini[i,1] == "XLLCORNER") { xn <- as.numeric(ini[i,2])
		  } else if (ini[i,1] == "XLLCENTER") { xc <- as.numeric(ini[i,2])
		  } else if (ini[i,1] == "YLLCORNER") { yn <- as.numeric(ini[i,2])
		  } else if (ini[i,1] == "YLLCENTER") { yc <- as.numeric(ini[i,2])
		  } else if (ini[i,1] == "CELLSIZE") { d <- as.numeric(ini[i,2])
		  } else if (ini[i,1] == "NODATA_VALUE") { try (nodataval <- as.numeric(ini[i,2]), silent=TRUE)
		  } else if (ini[i,1] == "NODATA") { try (nodataval <- as.numeric(ini[i,2]), silent=TRUE)		
		  }
    }  
	
	  if (is.na(nr)) stop('"NROWS" not detected') 
	  if (is.na(nc)) stop('"NCOLS" not detected')

	  if (is.na(nodataval)) {
	  	warning('"NODATA_VALUE" not detected. Setting it to -Inf\n  You can set it to another value with function "NAvalue"')
	  	nodataval <- -Inf
	  }
	
	  offwarn <- FALSE
	  if (is.na(d)) { 
	  	warning('"CELLSIZE" not detected. Setting it to 1.');
	  	offwarn = TRUE
	  	d <- 1 
	  } else if (d==0) {
	  	warning('"CELLSIZE" is reported as zero. Setting it to 1.');
	  	d <- 1 
	  }
	  d <- abs(d)

	  if (is.na(xn)) { 
	  	if (is.na(xc)) { 
	  		warning('"XLLCORNER" tag not detected. Setting it to 0.')
	  		offwarn = TRUE
	  		xn <- 0 
	  	} else {
	  		xn <- xn - 0.5 * d
	  	}
	  }
  	
  	if (is.na(yn)) { 
		  if (is.na(yc)) { 
		  	warning('"YLLCORNER" tag not detected. Setting it to 0.');
  			offwarn = TRUE
	  		yn <- 0
	  	} else {
	  		yn <- yc - 0.5 * d
	  	} 	
	  }
	  if (offwarn) {
	  	m <- 'The georeference of this object is probably wrong\n'
	  	if (offset != 6) {
	  		m <- paste(m, '  Are you using a wrong offset? Proceed with caution!\n', sep='')
	  	} 
	  	warning(m)
	  }
	
	  xx <- xn + nc * d
	  yx <- yn + nr * d
  } else {
    stop("File type not supported. Supported file types are: ESRI ASCII, GRASS ASCII")
  }

	x <- raster(ncols=nc, nrows=nr, xmn=xn, ymn=yn, xmx=xx, ymx=yx, crs='')	

	x@data@fromdisk <- TRUE
	x@file@offset <- offset
	x@file@driver <- 'ascii'
	x@file@nodatavalue <- nodataval
    x@file@name <- filename
	
	
	if (!is.na(crs)) {
		projection(x) <- crs
	}
	
  if(exists("multiplier")){
    x <- setValues(x, getValues(x)*multiplier)
  }
	return(x)
}


