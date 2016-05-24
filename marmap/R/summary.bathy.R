summary.bathy = function(object, ...){
	
	round(min(as.numeric(colnames(object))),2) -> lat.min
	round(max(as.numeric(colnames(object))),2) -> lat.max
	round(min(as.numeric(rownames(object))),2) -> lon.min
	round(max(as.numeric(rownames(object))),2) -> lon.max
	
	lon.max2 <- ifelse(lon.max > 180, lon.max-360, lon.max)
	flag.l1 <- ifelse(lon.min < 0, "W", "E")
	flag.l2 <- ifelse(lon.max2 < 0, "W", "E")

	flag.l3 <- ifelse(lat.min < 0, "S", "N")
	flag.l4 <- ifelse(lat.max < 0, "S", "N")

	one.minute = 0.016667
	as.numeric(rownames(object))[2] - as.numeric(rownames(object))[1] -> cell.size.centroid
	round(cell.size.centroid / one.minute, 1) -> cell.size.minute
	
	cat(paste("Bathymetric data of class 'bathy', with",dim(object)[1],"rows and",dim(object)[2],"columns\n"))
	cat(paste("Latitudinal range: ", lat.min," to ", lat.max, " (",abs(lat.min)," ",flag.l3," to ",abs(lat.max)," ",flag.l4,")\n",sep=""))
	cat(paste("Longitudinal range: ", lon.min," to ", lon.max, " (",abs(lon.min)," ",flag.l1," to ",abs(lon.max2)," ",flag.l2,")\n",sep=""))
	cat(paste("Cell size:",cell.size.minute,"minute(s)\n"))
	cat("\n")
	cat("Depth statistics:\n")
	print(summary(as.vector(object), ...))
	cat("\n")
	cat("First 5 columns and rows of the bathymetric matrix:\n")
	object[1:5, 1:5]
}