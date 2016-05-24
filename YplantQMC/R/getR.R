#'Get crown radius of a plant
#'
#'Calculates the crown radius of a plant, given the x,y,z coordinates of the
#'leaves. See Details. Typically not invoked by user, but used by
#'\code{\link{summary.plant3d}} to estimate crown width.
#'
#'The crown radius is the distance from the centre of the plant to the furthest
#'leaf in that direction. This distance is found for four directions
#'(corresponding to the four quadrants). The centre of the plant is found from
#'the average x,y coordinate (and may thus differ from the 'stem' location).
#'
#'@param xyz A matrix or dataframe of xyz coordinates.
#'@return Returns only the mean crown radius of the four quadrants.
#'@author Remko Duursma
#'@seealso
#'\code{\link{summary.plant3d}},\code{\link{crownhull}},\code{\link{Silhouette}}
#'@keywords misc
getR <- function(xyz){

	if(class(xyz) == "plant3d")
		xyz <- do.call("rbind", lapply(xyz$leaves, function(x)x$XYZ))
	
	# construct the hull (gives area and volume)
	if(is.list(xyz) && !is.data.frame(xyz))
		xyz <- as.matrix(do.call("rbind", xyz))
	else
		xyz <- as.matrix(xyz)

	# rescale
	x <- xyz[,1] - mean(xyz[,1])
	y <- xyz[,2] - mean(xyz[,2])
	z <- xyz[,3] - mean(xyz[,3])
	Xyz <- data.frame(x,y,z)
	r <- c()
	getr <- function(b){
		if(nrow(b) == 0)
			return(NA)
		else
			return(max(sqrt(b$x^2+b$y^2),na.rm=TRUE))
	}		
	r[1] <- getr(subset(Xyz, x>0 & y >0))
	r[2] <- getr(subset(Xyz, x>0 & y <0))
	r[3] <- getr(subset(Xyz, x<0 & y <0))
	r[4] <- getr(subset(Xyz, x<0 & y >0))
	
	return(mean(r,na.rm=TRUE))
}


# # To visualize: maybe should do this for all to check this algorithm!
# xyz <- XYZ[[21]]
# # rescale
# xyz$x <- xyz$x - mean(xyz$x)
# xyz$y <- xyz$y - mean(xyz$y)
# xyz$z <- xyz$z - mean(xyz$z)

# with(xyz, plot(x,y))
# r <- c()
# q1 <- subset(xyz, x>0 & y >0)
# r[1] <- with(q1, max(sqrt(x^2+y^2)))
# with(q1, points(x=c(0,x[which.max(sqrt(x^2+y^2))]),
# y=c(0,y[which.max(sqrt(x^2+y^2))]), type='l', col="red"))
# q1 <- subset(xyz, x>0 & y <0)
# r[2] <- with(q1, max(sqrt(x^2+y^2)))
# with(q1, points(x=c(0,x[which.max(sqrt(x^2+y^2))]),
# y=c(0,y[which.max(sqrt(x^2+y^2))]), type='l', col="red"))
# q1 <- subset(xyz, x<0 & y <0)
# r[3] <- with(q1, max(sqrt(x^2+y^2)))
# with(q1, points(x=c(0,x[which.max(sqrt(x^2+y^2))]),
# y=c(0,y[which.max(sqrt(x^2+y^2))]), type='l', col="red"))
# q1 <- subset(xyz, x<0 & y >0)
# r[4] <- with(q1, max(sqrt(x^2+y^2)))
# with(q1, points(x=c(0,x[which.max(sqrt(x^2+y^2))]),
# y=c(0,y[which.max(sqrt(x^2+y^2))]), type='l', col="red"))

