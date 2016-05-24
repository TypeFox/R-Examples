NCEP.loxodrome <- function(lat1,lat2,lon1,lon2){
# computes heading (in degrees from north) 
# between 2 positions on sphere
# Given: lat and lons in degrees
# J McLaren Dec 2008 ## translated to R by M. Kemp 14/09/2009
# see Gudmunsson, G A and T Alerstam ,1998, Optimal map projections for
# computing long-distance migration routes, J Avian Bio 29:597-605
# See also Alexander, J, 2004, Loxodromes: A Rhumb Way to Go,
# Mathematics Magazine 77/5:349-356 
# originally solved by Mercator graphically

## Multiplier for conversion from degrees to radians
deg2rad <- pi/180

## Create the Inverse cotangent function
acot <- function(x){
	return(atan(1/x))
}


# First transform to radians
lat1 <- deg2rad * lat1
lat2 <- deg2rad * lat2
lon1 <- deg2rad * lon1
lon2 <- deg2rad * lon2
deltaLon <- lon2 - lon1
 
# Then compute Mercator sigma co-ordinates
pi4 <- pi/4
Sig1 <- log(tan(pi4 + lat1/2))
Sig2 <- log(tan(pi4 + lat2/2))
deltaSig <- Sig2 - Sig1
 
# note if deltaSig < 0 latitude decreases ie bird travels south
#      if deltaLon < 0 longitude decreases ie bird travels west
 
# Now compute heading between locations 
# Need to know determine quadrant 
# matlab acot gives values from -90 to +90 (in radian equivalent)
 
if (deltaLon == 0 && deltaSig > 0){
    head <- 0 } else
if (deltaLon == 0 && deltaSig < 0){ 
    head <- 180 } else
if (deltaSig == 0 && deltaLon > 0){ 
    head <- 90 } else
if (deltaSig == 0 && deltaLon < 0){
    head <- 270 } else
if (deltaSig < 0 && deltaLon < 0){ 
    head <- acot(deltaSig/deltaLon)*180/pi + 180 } else
if (deltaSig < 0 && deltaLon > 0){
    head <- acot(deltaSig/deltaLon)*180/pi + 180 } else
if (deltaSig > 0 && deltaLon > 0){
    head <- acot(deltaSig/deltaLon)*180/pi } else
if (deltaSig > 0 && deltaLon < 0){
    head <- acot(deltaSig/deltaLon)*180/pi + 360 } else 
	stop("Migrant not moving, loxodrome = 'NA'") 

return (head)
}
