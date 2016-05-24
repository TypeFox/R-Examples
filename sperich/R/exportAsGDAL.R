exportAsGDAL <- 
function(grid, shift, resolution, directory=getwd(), filename="grid.tif", drivername="GTiff"){
	requireNamespace("rgdal")

	#check directory
	direc <- unlist(strsplit(directory,""))
	if (direc[length(direc)] != "/"){
		direc <- c(direc, "/")
	}
	directory <- ""
	for (i in 1:length(direc)){
		directory <- paste(directory,direc[i], sep="")
	}

	#get Dimension
	dimension <- dim(grid)

	#create dataframe
	values <- lat <- long <- vector(mode="numeric", length=dimension[1]*dimension[2])
	count <- 1
	for (m in 1:dimension[1]){
		for ( n in 1:dimension[2]){
			values[count] <- grid[m,n]
			long[count] <- shift[1] + resolution*(m-1)
			lat[count] <- shift[2] + resolution*(n-1)
			count <- count + 1
		}
	}
	result <- data.frame(values=values)
	points <- sp::SpatialPoints(data.frame(long=long, lat=lat), proj4string=sp::CRS("+proj=longlat +datum=WGS84"))
	#create SpatialGridDataFrame 
	result.spatial <- sp::SpatialPixelsDataFrame(points, result)	
	result.gdal <- as(result.spatial, "SpatialGridDataFrame")	

	rgdal::writeGDAL(result.gdal, fname=paste(directory,filename, sep=""), drivername=drivername)
}
