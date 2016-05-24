################################
#### Inverse of Lambert's equal aea projection
#### Tsagris Michail 3/2014 
#### mtsagris@yahoo.gr
################################

lambert.inv <- function(z, mu) {
  ## z contains the Lambert equal area projected data
  ## mu is the initial mean direction to which we will
  ## rotate the data after bringing them on the sphere
  z <- as.matrix(z)
  long <- ( atan(z[, 2] / z[, 1]) + pi * I(z[, 1] < 0) ) %% (2 * pi)
  lat <- 2 * asin( 0.5 * sqrt( rowSums(z^2) ) )
  u <- cbind(lat, long)  ## the data on the sphere in radians
  u <- u * 180 / pi  ## from radians to degrees
  y <- euclid(u)  ## the data in euclidean coordinates
  ## their mean direction is not exactly the north pole
  b <- c(0, 0, 1)  ## the north pole from which we will rotate the data
  mu <- mu / sqrt( sum(mu^2) )  ## make sure that mu is a unit vector
  H <- rotation(b, mu)  ## rotate the data so that their mean direction is mu
  u <- tcrossprod(y, H)
  euclid.inv(u)
}