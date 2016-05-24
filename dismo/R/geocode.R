# Author: Robert J. Hijmans
# License GPL3
# Version 1.0
# October 2010

geocode <- function(x, oneRecord=FALSE, extent=NULL, progress='', ...) {
	ntry <- list(...)$ntry
	if (is.null(ntry)) ntry <- 10
	reps <- min(ntry, 10)
	x <- as.character(x)
	xx <- unique(x)
	xx <- data.frame(ID=1:length(xx), place=xx)
	resall <- .geocode(xx$place, oneRecord=oneRecord, extent=extent, progress=progress)

	k <- which(is.na(resall[,3]))
	if (length(k) > 0 & reps > 1) {
		n <- 0
		for (i in 2:reps) {
			print(paste('try', i, '...'))
			flush.console()
			j <- which(is.na(resall[,2]))
			if (length(j) == 0) break
			ids <- unique(resall[j,1])
			res <- .geocode(xx[ids, 'place'], oneRecord=oneRecord, extent=extent, progress=progress)
			k <- which(!is.na(res[,2]))
			if (length(k) == 0) {
				if (n == 2) break
				n <- n + 1
			} else {
				n <- 0
			}
			if (oneRecord) {
				res$ID <- j[match(res$ID, 1:length(j))]
				resall[j,] <- res
			} else {
				resall <- resall[-j, ,drop=FALSE]
				res$ID <- j[match(res$ID, 1:length(j))]
				resall <- rbind(resall, res)
			}
		}
	}
	x <- data.frame(ID=1:length(x), originalPlace=x)
	resall$ID <- NULL
	res <- merge(x, resall, by='originalPlace', all.x=TRUE)
	rownames(res) <- NULL
	res <- res[order(res$ID), , drop=FALSE]
	res$ID <- NULL
	res
}


.geocode <- function(x, oneRecord=FALSE, extent=NULL, progress='', ...) {
	
	if (! requireNamespace('XML')) stop('You need to install the XML package to be able use this function')

	burl <- "http://maps.google.com/maps/api/geocode/xml?address="
	res1 <- data.frame(matrix(NA, ncol=8, nrow=1))
	colnames(res1) <- c('ID', 'interpretedPlace', 'longitude', 'latitude', 'xmin', 'xmax', 'ymin', 'ymax')
	res <- res1[-1, ,drop=FALSE]
	pb <- pbCreate(length(x), progress)
	for (z in 1:length(x)) {
		r <- x[z]
		r <- gsub(', ', ',', r)
		r <- gsub(' ,', ',', r)
		r <- trim(r)
		if (length(r) > 0 & !is.na(r)) {
			r <- gsub(' ', '+', r)
			if (is.null(extent)) {
				gurl <- paste(burl, r, "&sensor=false", sep="")
			} else {
				e <- extent(extent)
				extent <- paste(e@ymin,',',e@xmin,'|',e@ymax,',',e@xmax,sep='')
				gurl <- paste(burl, r, "&bounds=", extent, "&sensor=false", sep="")	
				gurl <- iconv(gurl, to='UTF-8')
			}
			try( doc <- XML::xmlInternalTreeParse(gurl, isURL=TRUE) )
			if (class(doc)[1] == 'try-error') {
				warning('cannot parse XML document\n')
				status <- ''
			} else { 
				status <- XML::xmlValue(XML::getNodeSet(doc, "//GeocodeResponse//status")[[1]])
			}
			
			if (status != "OK") {
				message(status, ':', r)
				w <- res1 
				w[1] <- z
				res <- rbind(res, w)
				next
			}
		
			p <- XML::xmlToList(doc)
			n <- length(p)-1
			place <- rep(NA, n)
			location <- matrix(ncol=2, nrow=n)
			viewport <- matrix(ncol=4, nrow=n)
			bounds <- matrix(ncol=4, nrow=n)
			for (i in 1:n) {
				plc <- p[i+1]$result$formatted_address
				place[i] <- ifelse(is.null(plc), "", plc)
				location[i,] <- as.numeric(c(p[i+1]$result$geometry$location$lng, p[i+1]$result$geometry$location$lat))
				viewport[i,] <- as.numeric(c(p[i+1]$result$geometry$viewport$southwest$lng, p[i+1]$result$geometry$viewport$northeast$lng, p[i+1]$result$geometry$viewport$southwest$lat, p[i+1]$result$geometry$viewport$northeast$lat) )
				bnds <- as.numeric(c(p[i+1]$result$geometry$bounds$southwest$lng, p[i+1]$result$geometry$bounds$northeast$lng, p[i+1]$result$geometry$bounds$southwest$lat, p[i+1]$result$geometry$bounds$northeast$lat) )
				if (length(bnds)==4) bounds[i,] <- bnds
			}

			w <- cbind(viewport, bounds)
			w[,c(1,3)] <- pmin(w[,c(1,3)], w[,c(1,3)+4], na.rm=TRUE)
			w[,c(2,4)] <- pmax(w[,c(2,4)], w[,c(2,4)+4], na.rm=TRUE)
			w <- w[,1:4,drop=FALSE]
			if (oneRecord & nrow(w) > 1) {
				f <- apply(w, 2, range)
				g <- apply(location, 2, mean)
				w <- data.frame(z, NA, g[1], g[2], f[1,1], f[2,2], f[1,3], f[2,4])
			} else {
				w <- data.frame(z, place, location, w)
			}
			colnames(w) <- colnames(res)
			res <- rbind(res, w)
		} else {
			w <- res1
			w[1] <- z
			colnames(w) <- colnames(res)
			res <- rbind(res, w)
		}
		pbStep(pb, z) 
	} 
	pbClose(pb)

	rownames(res) <- 1:nrow(res)
	
	da <- pointDistance(res[,3:4], res[,c(5,7)], longlat=T)
	db <- pointDistance(res[,3:4], res[,c(6,8)], longlat=T)
	res$uncertainty <- round(pmin(da, db))
	xx <- data.frame(ID=1:length(x), originalPlace=x)
	
	merge(res, xx, by='ID')
}


#a = geocode('San Jose, Mexico', oneRecord=F)

 
 