################################
#### Change between latitude/longitude and Cartesian coordinates
#### Tsagris Michail 1/2014 
#### mtsagris@yahoo.gr
################################

euclid.inv <- function(U) {
  ## U is a 3-column matrix of unit vectors
  ## the cartesian coordinates
  U <- as.matrix(U)
  if (ncol(U) == 1)  U <- t(U)
  u <- cbind( acos(U[, 1]), ( atan(U[, 3]/U[, 2]) + pi * I(U[, 2]<0) )
  %% (2 * pi) )
  u <- u * 180/pi  ## from rads to degrees
  colnames(u) <- c("Lat", "Long")
  ## u is a matrix of two columns
  ## the first column is the latitude and the second the longitude in degrees
  u
}