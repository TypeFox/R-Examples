################################
#### Change between latitude/longitude and Cartesian coordinates
#### Tsagris Michail 1/2014 
#### mtsagris@yahoo.gr
################################

euclid <- function(u) {
  ## u is a matrix of two columns
  ## the first column is the latitude and the second the longitude
  u <- as.matrix(u)
  if (ncol(u) == 1)  u <- t(u)
  u <- pi * u/180  ## from degrees to rads
  U <- cbind(cos(u[, 1]), sin(u[, 1]) * cos(u[, 2]), sin(u[, 1]) * sin(u[, 2]))
  colnames(U) <- c("x", "y", "z")
  ## U are the cartesian coordinates of u
  U
}