#'Calculates the area of the 2D convex hull of a projected plant
#'
#'@description Calculates the so-called 'silhouette area' of a projected plant. That is, the area of the (convex) outline of a projection of a plant onto a plane.
#'
#'
#'@param obj Either a \code{plant3d} object, or a \code{projectedplant3d}
#'object.
#'@param azimuth,altitude Viewing azimuth and altitude (ignored if a projected
#'plant is given as input).
#'@note Not usually called by the user. Use instead, \code{\link{STARbar}} and
#'\code{\link{projectplant}}. For the latter, see the option
#'\code{silhouette=TRUE}.
#'@author Remko Duursma
#'@seealso
#'\code{\link{summary.plant3d}},\code{\link{STARbar}},\code{\link{projectplant}}
#'@keywords misc
#'@examples
#'
#'
#'# Silhouette returns the area of the 2D convex hull (H), and the coordinates of the 2D hull (xyz):
#'Silhouette(sugarmaple, altitude=0, azimuth=45)
#'
#'@export
Silhouette <- function(obj, azimuth=NA, altitude=NA){

	if(inherits(obj, "plant3d")){
		if(is.na(azimuth) | is.na(altitude))stop("Need azimuth and altitude viewing angles.")
		newplant <- projectplant(obj, azimuth, altitude)
     } else	
		newplant <- obj
		
	P <- do.call("rbind", lapply(newplant$leaves, function(x)x$XYZ[,1:2]))
	ch <- chull(P)
	chp <- c(ch, ch[1])
	H <- areapoly(P[chp,])
	
	l <- list()
	l$H <- H
	l$xyz <- P[chp,]
	
return(l)
}

