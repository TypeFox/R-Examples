Kn <- function(plant, kneighbors=5, returnmatrix=FALSE){

	# Coordinates
	if(class(plant) == "plant3d"){
		Xyz <- list(plant$leafbasecoor, plant$leaftipcoor)
		
		# Get coordinates of leaf midpoints
		xyz <- data.frame(x = (Xyz[[1]][,1] + Xyz[[2]][,1])/2,
					  y = (Xyz[[1]][,2] + Xyz[[2]][,2])/2,
					  z = (Xyz[[1]][,3] + Xyz[[2]][,3])/2 )
		nleaves <- plant$nleaves
		
	} else {
		Xyz <- plant
		xyz <- data.frame(x=Xyz[,1], y=Xyz[,2], z=Xyz[,3])
		nleaves <- nrow(xyz)
	}
		
  x <- xyz$x / 1000
	y <- xyz$y / 1000
	z <- xyz$z / 1000
	
	res <- matrix(nrow=nleaves,ncol=nleaves)
	res[] <- 0
	
	f <- .Fortran("kn", as.double(x),as.double(y),as.double(z),as.integer(nleaves),as.double(res),PACKAGE="YplantQMC")

	m <- f[[5]]
	m <- matrix(m, ncol=nleaves)
	
	if(returnmatrix)return(m)
	
	m[m == 0] <- NA
	
	meanminn <- function(x,k)mean(sort(x)[1:k], na.rm=TRUE)
	if(length(kneighbors) == 1){
		o <- apply(m, 1, meanminn,kneighbors)
		return(mean(o))
	} else {
		o <- list()
		for(j in 1:length(kneighbors)){
			if(j > (nleaves-1)){o[[j]] <- NA; next}
			o[[j]] <- apply(m, 1, meanminn,kneighbors[j])
		}
		return(sapply(o,mean,na.rm=TRUE))
	}
}

