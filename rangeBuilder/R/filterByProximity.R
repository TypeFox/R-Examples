# Function to filter out points that are closer than some distance from each other

filterByProximity <- function(xy, dist, mapUnits = FALSE, returnIndex = FALSE) {

	#xy can be either a SpatialPoints or SPDF object, or a matrix
	#dist is in km if mapUnits=F, in mapUnits otherwise
	#returnIndex = TRUE will return indices of points that would be dropped
	if (class(xy) == 'data.frame') {
		xy <- as.matrix(xy)
	}
	if (!mapUnits) {
		d <- spDists(xy,longlat = TRUE)
	}
	if (mapUnits) {
		if (is.matrix(xy)) {
			d <- as.matrix(dist(xy, method = 'euclidian'))
		} else {
			d <- spDists(xy, longlat = FALSE)
		}
	}
	diag(d) <- NA
	close <- d <= dist
	diag(close) <- NA
	closePts <- which(close, arr.ind = TRUE)
	discard <- matrix(nrow = 2, ncol=2)
	if (nrow(closePts) > 0) {
		while (nrow(closePts) > 0) {
			if ((!paste(closePts[1,1],closePts[1,2],sep='_') %in% paste(discard[,1],discard[,2],sep='_')) & (!paste(closePts[1,2],closePts[1,1],sep='_') %in% paste(discard[,1],discard[,2],sep='_'))) {
				discard <- rbind(discard, closePts[1,])
				closePts <- closePts[-union(which(closePts[,1] == closePts[1,1]), which(closePts[,2] == closePts[1,1])),]
			}
		}
		discard <- discard[complete.cases(discard),]
		if (class(discard) == 'matrix') {
			if (returnIndex) {
				return(discard[,1])
			} else {
				return(xy[-discard[,1],])
			}
		} else {
			if (returnIndex) {
				return(discard[1])
			} else {
				return(xy[-discard[1],])
			}
		}
	}
	if (nrow(closePts) == 0) {
		if (returnIndex) {
			return(NA)
		} else {
			return(xy)
		}
	}
}

