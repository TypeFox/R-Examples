path.profile <- function(path,bathy,plot=FALSE,...) {
	
	if (class(bathy)!="bathy") stop("Object bathy is not of class bathy")
	if (ncol(path)!=2) stop("Object path should have 2 columns: Longitude and Latitude")
	
	out <- matrix(0,ncol=4)
	colnames(out) <- c("lon","lat","dist.km","depth")
	
	for (i in 1:(nrow(path)-1)) {
		df <- get.transect(mat=bathy, x1=path[i,1], y1=path[i,2],x2=path[i+1,1],y2=path[i+1,2],distance=TRUE)
		# df <- df[-1,]
		df[,3] <- df[,3] + out[nrow(out),3]
		out <- rbind(out,df)
	}

	out <- unique(out[-1,])

	if (plot){
		dev.new()
		plotProfile(out,...)
	} 
	return(out)
}

