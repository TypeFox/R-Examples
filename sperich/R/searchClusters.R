searchClusters <-
function(species.richness, dimension, shift, resolution, clusterlimit){
	requireNamespace("SDMTools")
	requireNamespace("sp")

	#create a binary matrix out of species.richness
	grid <- matrix(0, dimension[1], dimension[2])
	#fill binary matrix
	grid[which(species.richness < clusterlimit)] <- 0
	grid[which(species.richness >= clusterlimit)] <- 1

	#SpatialGridDataFrame needed, first step: create dataframe 
	values <- vector(mode="integer", length=dimension[1]*dimension[2])
	lat <- long <- vector(mode="numeric", length=dimension[1]*dimension[2])
	count <- 1
	for (m in 1:dimension[1]){
		for ( n in 1:dimension[2]){
			values[count] <- grid[m,n]
			long[count] <- shift[1] + resolution*n
			lat[count] <- shift[2] + resolution*m
			count <- count + 1
		}
	}
	result <- data.frame(values=values)
	points <- sp::SpatialPoints(data.frame(long=long, lat=lat), proj4string=sp::CRS("+proj=longlat +datum=WGS84"))
	#create SpatialGridDataFrame 
	result.spdf <- sp::SpatialPixelsDataFrame(points, result)
	result.sgdf <- as(result.spdf, "SpatialGridDataFrame")	

	#start clustering
	cluster <- SDMTools::ConnCompLabel(result.sgdf)

	#convert dataframe to matrix
	clustermatrix <- as.matrix(cluster)
	#adjust matrix
	clustermatrix <- t(clustermatrix)[dim(t(clustermatrix))[1]:1,]	

	#create list for clusters
	clusterlist <- list()	

	#fill list with clusters
	if (max(clustermatrix)==0){
		clusterlist[[1]] <- which(clustermatrix==0)
	} else {
		for (i in 1:max(clustermatrix)){
			clusterlist[[length(clusterlist)+1]] <- which(clustermatrix==i)
		}
	}
	return(clusterlist)
}
