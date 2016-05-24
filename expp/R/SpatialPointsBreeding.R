
SpatialPointsBreeding <- function(data, 
                                  proj4string, 
                                  coords = ~ x + y, 
                                  breeding = ~ male + female, 
                                  id ) {
	d = data
	row.names(d) = NULL
	d$k = 1:nrow(d)
	coordinates(d) <- coords
  if(missing(proj4string)) proj4string  = CRS(as.character(NA))
	proj4string(d) = proj4string
	
	ids = data[, id]
	
	m = as.character(breeding[[2]][2])
	f  = as.character(breeding[[2]][3])
	males = data[, m]
	females = data[, f]
	
	d@data[, m]  = NULL
	d@data[, f]  = NULL
	d@data[, id] = NULL


	
	new("SpatialPointsBreeding", d, id = ids, male = males, female= females)
}



if (!isGeneric("plot")) setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

	
setMethod("plot", signature(x = "SpatialPointsBreeding", y = "missing"),
          function(x, pch = 20, axes = FALSE, add = FALSE, 
                   xlim = NULL, ylim = NULL, ..., cex = 1, col = "grey", lwd = 1, bg = "grey90") {
            if (! add)
              plot(as(x, "Spatial"), axes = axes, xlim = xlim, ylim = ylim, ...)
            cc = coordinates(x)
            points(cc[,1], cc[,2], pch = pch, cex = cex, col = col, lwd = lwd, bg = bg)
            text(cc[,1], cc[,2], x@id, pos = 4, cex = cex)
            text(cc[,1], cc[,2], x@female, pos = 1,  cex =  cex-0.1)
            text(cc[,1], cc[,2], x@male, pos = 3,  cex = cex-0.1)
            
          })
		  

setMethod("plot", signature(x = "SpatialPointsBreeding", y = "eppMatrix"),
          function(x, y, pch = 20, axes = FALSE, add = FALSE, 
                   xlim = NULL, ylim = NULL, ..., cex = 1, col = "grey", col.epp = "red", lwd = 1, lty = 2, 
                   bg = "grey90") {
            if (! add)
              plot(as(x, "Spatial"), axes = axes, xlim = xlim, ylim = ylim, ...)
            cc = coordinates(x)
			
			# nests
			points(cc[,1], cc[,2], pch = pch, cex = cex, col = col, lwd = lwd, bg = bg)
            text(cc[,1], cc[,2], x@id, pos = 4, cex = cex)
            
			# males
			epm = which(x@male %in% y@male)
			try(text(cc[epm,1], cc[epm,2], x@male[epm], pos = 1,  cex =  cex-0.1, col = col.epp), silent = TRUE)
			try(text(cc[-epm,1], cc[-epm,2], x@male[-epm], pos = 1,  cex =  cex-0.1), silent = TRUE)
			
			# females
			epf = which(x@female %in% y@female)
			try(text(cc[epf,1], cc[epf,2], x@female[epf], pos = 3,  cex =  cex-0.1, col = col.epp), silent = TRUE)
			try(text(cc[-epf,1], cc[-epf,2], x@female[-epf], pos = 3,  cex =  cex-0.1), silent = TRUE)
			
			# connections
			for(i in 1:length(y@male) ) {
				mc = cc[which(x@male == y@male[i]), ]
				if(length(mc) == 0) warning("EP male ", sQuote(y@male[i]), " not found.")
				if( length(mc) > 2) mc = mc[1, ] # only one line per polygynous male
					
				fc = cc[which(x@female == y@female[i]), ]
				if(length(fc) == 0) warning("EP female ", sQuote(y@female[i]), " not found.")
								
				arrows(mc[1], mc[2], fc[1], fc[2], col = col.epp, code = 3, angle = 12, length = 0.2, lty = lty)
				}

			
			
			
			
            
          })
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  

	











