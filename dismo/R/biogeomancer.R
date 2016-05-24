# Download geographic data and return as R object
# Author: Robert J. Hijmans, r.hijmans@gmail.com
# License GPL3
# Version 0.1
# October 2008


.biogeomancer <- function(country='', adm1='', adm2='', locality='', singleRecord=TRUE, progress='text') {

#	if (! require(XML)) stop('You need to install the XML package to be able use this function')

	d <- data.frame(country, adm1, adm2, locality)
	d[is.na(d)] <- ''
	pb <- pbCreate(dim(d)[1], progress)

	for (r in 1:dim(d)[1]) {
		theurl <- paste("http://bg.berkeley.edu:8080/ws/single?cy=", d$country[r], "&sp=", d$adm1[r], "&co=", d$adm2[r], "&locality=", d$locality[r], sep='')
		
		try( doc <- XML::xmlInternalTreeParse(theurl) )
		if (class(doc)[1] == 'try-error') {
			ans <- data.frame(lon=NA, lat=NA, coordUncertaintyM=NA)
		} else {
# to do: improved parsing:	
			nodes <- XML::getNodeSet(doc, "//georeference")
			if(length(nodes) == 0) {
				ans <- data.frame(lon=NA, lat=NA, coordUncertaintyM=NA)
			} else {
				varNames <- c("decimalLongitude", "decimalLatitude", "geodeticDatum", "coordinateUncertaintyInMeters")
				dims <- c(length(nodes), length(varNames)) 
   # create an empty data frame with as many rows and columns as needed.
				ans <- as.data.frame(replicate(dims[2], rep(as.character(NA), dims[1]), simplify = FALSE), stringsAsFactors = FALSE)
				names(ans) <- varNames
    # Fill in the rows based on the names.
				for(i in seq(length = dims[1])) { 
					ans[i, varNames] = XML::xmlSApply(nodes[[i]], XML::xmlValue)[varNames]
				}
				ans <- ans[,-3]
				names(ans) <- c("lon", "lat", "coordUncertaintyM")
				if (singleRecord) {
					ans <- ans[which.min(ans[,"coordUncertaintyM"]),][1,]
				}
			}
		}
		if (r == 1) {
			res <- cbind(id=r, ans)
		} else {
			res <- rbind(res, cbind(id=r, ans))
		}
		pbStep(pb, r) 
	} 
	pbClose(pb)
	res$lon = as.numeric(res$lon)
	res$lat = as.numeric(res$lat)
	res$coordUncertaintyM = as.numeric(res$coordUncertaintyM)
	res
}

# biogeomancer(country='United States', adm1='California', adm2='Yolo', locality=c('Davis', 'Woodland'))

