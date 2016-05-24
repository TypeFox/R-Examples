# Author: Robert Hijmans
# November 2009
# License GPL3
#
# based on:
# MkBCvars.AML 
# Author Robert Hijmans
# January 2006  
# Museum of Vertebrate Zoology, UC Berkeley
#
# Version 2.3
#
# function to create19 BIOCLIM variables from 
# monthly T-min, T-max, and Precipitation data
#
# BIO1 = Annual Mean Temperature
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# BIO3 = Isothermality (P2/P7) (* 100)
# BIO4 = Temperature Seasonality (standard deviation *100)
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO7 = Temperature Annual Range (P5-P6)
# BIO8 = Mean Temperature of Wettest Quarter
# BIO9 = Mean Temperature of Driest Quarter
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO12 = Annual Precipitation
# BIO13 = Precipitation of Wettest Month
# BIO14 = Precipitation of Driest Month
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
# BIO16 = Precipitation of Wettest Quarter
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter
# 
# These summary Bioclimatic variables are after:
#   Nix, 1986. A biogeographic analysis of Australian elapid snakes. In: R. Longmore (ed.).
#      Atlas of elapid snakes of Australia. Australian Flora and Fauna Series 7.
#      Australian Government Publishing Service, Canberra.
#
# and Expanded following the ANUCLIM manual
# 


if (!isGeneric("biovars")) {
	setGeneric("biovars", function(prec, tmin, tmax, ...)
		standardGeneric("biovars"))
}	


setMethod('biovars', signature(prec='vector', tmin='vector', tmax='vector'), 
	function(prec, tmin, tmax) {
		biovars(t(as.matrix(prec)), t(as.matrix(tmin)), t(as.matrix(tmax)))
	}
)


setMethod('biovars', signature(prec='Raster', tmin='Raster', tmax='Raster'), 
	function(prec, tmin, tmax, filename='', progress='', ...) {

	if (nlayers(prec) != 12) stop('nlayers(prec) is not 12')
	if (nlayers(tmin) != 12) stop('nlayers(tmin) is not 12')
	if (nlayers(tmax) != 12) stop('nlayers(tmax) is not 12')
	
	
	# temporary fix to avoid warning
	compareRaster(prec, tmin, tmax)

	out <- brick(prec, values=FALSE)
	out@data@nlayers <- as.integer(19)
	names(out) <- paste('bio', 1:19, sep="")
	
	filename <- trim(filename)
	if (!canProcessInMemory(out, 18)) {
		if (filename == '') { 
			filename <- rasterTmpFile()
		}
	} 
	if (filename == "") {
		v <- matrix(nrow=ncell(out), ncol=19)
	} else {
		out <- writeStart(out, filename, ...)
	}	

	tr <- blockSize(out, n=nlayers(out)+36)
	pb <- pbCreate(tr$n, ...)	
	for (i in 1:tr$n) {
		prc <- getValues(prec, tr$row[i], tr$nrows[i])
		tmn <- getValues(tmin, tr$row[i], tr$nrows[i])
		tmx <- getValues(tmax, tr$row[i], tr$nrows[i])
		p <- biovars(prc, tmn, tmx)
		if (filename != "") {
			out <- writeValues(out, p, tr$row[i])
		} else {
			start <- (tr$row[i]-1) * out@ncols + 1
			end <- (tr$row[i]+tr$nrows[i]-1) * out@ncols
			v[start:end,] <- p
		}
	}
	
	if (filename == "") {
		out <- setValues(out, v)
	} else {
		out <- writeStop(out)
	}	
	return(out)
}
)




setMethod('biovars', signature(prec='matrix', tmin='matrix', tmax='matrix'), 
	function(prec, tmin, tmax) {

		if (nrow(prec) != nrow(tmin) | nrow(tmin) != nrow(tmax) ) {
			stop('prec, tmin and tmax should have same length')
		}
		
		if (ncol(prec) != ncol(tmin) | ncol(tmin) != ncol(tmax) ) {
			stop('prec, tmin and tmax should have same number of variables (columns)')
		}
		
		# can't have missing values in a row
		nas <- apply(prec, 1, function(x){ any(is.na(x)) } )
		nas <- nas | apply(tmin, 1, function(x){ any(is.na(x)) } )
		nas <- nas | apply(tmax, 1, function(x){ any(is.na(x)) } )
		p <- matrix(nrow=nrow(prec), ncol=19)
		colnames(p) = paste('bio', 1:19, sep='')
		if (all(nas)) { return(p) }
		
		prec[nas,] <- NA
		tmin[nas,] <- NA
		tmax[nas,] <- NA
		
		window <- function(x)  { 
			lng <- length(x)
			x <- c(x,  x[1:3])
			m <- matrix(ncol=3, nrow=lng)
			for (i in 1:3) { m[,i] <- x[i:(lng+i-1)] }
			apply(m, MARGIN=1, FUN=sum)
		}
		
		tavg <- (tmin + tmax) / 2
		
# P1. Annual Mean Temperature 
		p[,1] <- apply(tavg,1,mean)
# P2. Mean Diurnal Range(Mean(period max-min)) 
		p[,2] <- apply(tmax-tmin, 1, mean)
# P4. Temperature Seasonality (standard deviation) 
		p[,4] <- 100 * apply(tavg, 1, sd)
# P5. Max Temperature of Warmest Period 
		p[,5] <- apply(tmax,1, max)
# P6. Min Temperature of Coldest Period 
		p[,6] <- apply(tmin, 1, min)
# P7. Temperature Annual Range (P5-P6) 
		p[,7] <- p[,5] - p[,6]
# P3. Isothermality (P2 / P7) 
		p[,3] <- 100 * p[,2] / p[,7]
# P12. Annual Precipitation 
		p[,12] <- apply(prec, 1, sum)
# P13. Precipitation of Wettest Period 
		p[,13] <-  apply(prec, 1, max)
# P14. Precipitation of Driest Period 
		p[,14] <-  apply(prec, 1, min)
# P15. Precipitation Seasonality(Coefficient of Variation) 
# the "1 +" is to avoid strange CVs for areas where mean rainfaill is < 1)
		p[,15] <- apply(prec+1, 1, cv)
		
# precip by quarter (3 months)		
		wet <- t(apply(prec, 1, window))
# P16. Precipitation of Wettest Quarter 
		p[,16] <- apply(wet, 1, max)
# P17. Precipitation of Driest Quarter 
		p[,17] <- apply(wet, 1, min)
		tmp <- t(apply(tavg, 1, window)) / 3
		
		if (all(is.na(wet))) {
			p[,8] <- NA		
			p[,9] <- NA		
		} else {
# P8. Mean Temperature of Wettest Quarter 
			wetqrt <- cbind(1:nrow(p), as.integer(apply(wet, 1, which.max)))
			p[,8] <- tmp[wetqrt]
# P9. Mean Temperature of Driest Quarter 
			dryqrt <- cbind(1:nrow(p), as.integer(apply(wet, 1, which.min)))
			p[,9] <- tmp[dryqrt]
		}
# P10 Mean Temperature of Warmest Quarter 
		p[,10] <- apply(tmp, 1, max)

# P11 Mean Temperature of Coldest Quarter
		p[,11] <- apply(tmp, 1, min) 

		if (all(is.na(tmp))) {
			p[,18] <- NA		
			p[,19] <- NA
		} else {
# P18. Precipitation of Warmest Quarter 
			hot <- cbind(1:nrow(p), as.integer(apply(tmp, 1, which.max)))
			p[,18] <- wet[hot]
# P19. Precipitation of Coldest Quarter 
			cold <- cbind(1:nrow(p), as.integer(apply(tmp, 1, which.min)))
			p[,19] <- wet[cold]
		}
		
		return(p)	
	}
)

