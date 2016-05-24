DAPA.plant3d <- function (object, azimuth=NA, altitude=NA, exact=FALSE, npside=50, 
                      returnpolygons=FALSE, quiet=FALSE, progressbar=!quiet,
					  ...) 
{

	plant <- object
	AZ <- azimuth * 180/pi
	ALT <- altitude * 180/pi
	
  if(exact)
    stop("Method 'exact' temporarily unavailable. Contact package maintainer.")

	# Get total leaf area.
	LA <- plant$leafarea * 10^6
	nleaves <- plant$nleaves
	
	# Project plant
  newplant <- projectplant(plant, AZ, ALT)
    
	# Silhouette area (area of convex hull in 2D projected plane).
	silho <- Silhouette(newplant)$H
	
	# List with XY vertices of leaves on projected plane.
	P <- lapply(newplant$leaves, function(x)x$XYZ[,1:2])
		
	# Projected area (two methods give same answer:
	# The first is the area of the projected polygon, the second uses the acosangle
	# between leaf normal and viewing direction (a bit faster).
	PA <- sum(sapply(P, areapoly))
	#PA <- sum(newplant$acosangle * plant$leafdata$area)

# 	if (exact) {
# 		allunion <- fastunion(P)
# 		DA <- area.poly(allunion)
# 	}
# 	if (!exact) {
		ray <- gridtrace(newplant, npside = npside, returnall=FALSE)
		DA <- ray$nintersect * ray$pointarea
# 	}
	
	l <- list(DA = DA, PA = PA, azimuth = AZ, altitude = ALT, 
			exact = exact, LA = LA, H=silho)

return(l)
}
