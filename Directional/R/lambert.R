################################
#### Lambert's equal area projection
#### Tsagris Michail 3/2014 
#### mtsagris@yahoo.gr
#### References: Kent John T. (1982, JRSS Series C)
#### The Fisher-Bingham distribution on the sphere. 
################################

lambert <- function(y) {
  ## y contains the data in degrees, latitude and longitude
  u <- euclid(y)  ## transform the data into euclidean coordinates
  m <- colMeans(u)
  m <- m / sqrt(sum( m^2) )  ## the mean direction
  b <- c(0, 0, 1)  ## the north pole 
  H <- rotation(m, b)  ## the rotation matrix
  u1 <- tcrossprod(u, H)  ## rotating the data so that their mean 
  ## direction is the north pole
  u2 <- euclid.inv(u1)  ## bring the data into degrees again
  u2 <- pi * u2 / 180  ## from degrees to radians
  theta <- u2[, 1]
  phi <- u2[, 2]
  rho <- 2 * sin(theta / 2)  ## radius of the disk is sqrt(2)
  z1 <- rho * cos(phi)  ## x coordinate
  z2 <- rho * sin(phi)  ## y coordinate
  cbind(z1, z2)  ## the Lambert equal area projected data on the disk
}