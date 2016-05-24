################################################################################
## Computes the distance along the orthodrome in km between two points on the Earth
## x1, x2: longitudes of points 1 and 2 in decimal coordinates, W = negative sign
## y1, y2: latitudes  of points 1 and 2 in decimal coordinates, S = negative sign
## Reference: Imboden, C. & D. Imboden (1972). Vogelwarte 26: 336-346.
## Author: F. Korner, Jan 2014
## depends on package geosphere
#################################################################################
orthodrom.dist<-function(x1, y1, x2, y2){
  dist <- distMeeus(cbind(x1, y1), cbind(x2, y2))/1000
  return(dist)
}
###################################################################################