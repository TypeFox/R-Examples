zip.plot <- function(data, zip.file = system.file("data", "zips.tab", package = "muRL"), map.type = "state", cex = 1, col = "black", pch = 20, jitter.factor = NULL, ...){
	
	data$zip <- as.numeric(substr(data$zip, 1, 5))
	
	zips <- read.delim(zip.file, stringsAsFactor = FALSE)
	
	data.z <- merge(data, zips[,c("zip", "lat", "lon")], by.x = "zip", by.y = "zip", all.x = TRUE)

	maps::map(map.type, ...)
	if(!(is.null(jitter.factor))){
		data.z$lon <- jitter(data.z$lon, factor = jitter.factor)
		data.z$lat <- jitter(data.z$lat, factor = jitter.factor)
	}
	points(data.z$lon, data.z$lat, cex = cex, col = col, pch = pch)
}