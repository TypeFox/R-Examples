# function to change the sign of long and lat and check against country
# currently implemented for longlat WGS84 (as this is the coord system of the world points)

flipSign <- function(coordVec, country, returnMultiple = FALSE, proj = "+proj=longlat +datum=WGS84") {
#coordVec is a vector: long, lat
# country is the name of the country associated with those coordinates
# returnMultiple: if multiple sign flips lead to the correct country, return all options. If FALSE, returns the coords with the fewest needed sign flips.

	# check that provided proj4string is valid
	CRS(proj)

	if (is.matrix(coordVec) | is.data.frame(coordVec)) {
		coordVec <- as.numeric(coordVec)
	}

	country <- toupper(country)

	if (!country %in% unlist(countryList)) {
		stop('Country is not recognized.')
	}

	if (country %in% closestCountry(coordVec, proj = proj)) {
		cat('\tCoordinates already match country.\n')
		names(coordVec) <- c('long','lat')
		return(list(matched = TRUE, newcoords = coordVec))
	}

	long <- coordVec[1]
	lat <- coordVec[2]

	# create alternative signs
	allcoords <- matrix(ncol = 2, nrow = 7)
	colnames(allcoords) <- c('long','lat')
	allcoords[1,] <- c(long*-1, lat)
	allcoords[2,] <- c(long, lat*-1)
	allcoords[3,] <- c(lat, long)
	allcoords[4,] <- c(lat*-1, long)
	allcoords[5,] <- c(lat, long*-1)
	allcoords[6,] <- c(long*-1, lat*-1)
	allcoords[7,] <- c(lat*-1, long*-1)
	
	#but remove nonsensical coordinates
	bb <- SpatialPoints(matrix(c(-180, -90, 180, 90), nrow=2, ncol=2, byrow=TRUE), proj4string=CRS('+proj=longlat +datum=WGS84'))
	bb <- sp::spTransform(bb, CRS(proj))
	bb <- bbox(bb)

	drop <- union(which(allcoords[,'lat'] > bb[2,2]), which(allcoords[,'lat'] < bb[2,1]))
	if (length(drop) > 0) {allcoords <- allcoords[-drop,]}

	allcountries <- apply(allcoords, 1, function(x) closestCountry(x, proj = proj))

	if (country %in% unlist(allcountries)) {
		match <- lapply(allcountries, function(x) country %in% x)
		newcoords <- allcoords[which(match == TRUE),]
		newcoords <- matrix(newcoords, ncol = 2, byrow = TRUE)
		colnames(newcoords) <- c('long','lat')

		if (nrow(newcoords) == 2 & returnMultiple == FALSE) {
			if (identical(as.numeric(c(newcoords[1,1] * (-1), newcoords[1,2])), coordVec) | identical(as.numeric(c(newcoords[1,1], newcoords[1,2] * (-1))), coordVec)) {
				newcoords <- newcoords[1,]
			} else if (identical(as.numeric(c(newcoords[2,1] * (-1), newcoords[2,2])), coordVec) | identical(as.numeric(c(newcoords[2,1], newcoords[2,2] * (-1))), coordVec)) {
				newcoords <- newcoords[2,]
			}
		}

		return(list(matched = TRUE, newcoords = newcoords))
	} else {
		return(list(matched = FALSE, newcoords = NA))
	}
}

# if two coordinates return the correct country, we will return the one with the fewest changes needed