################################
#### Rotation axis and angle of rotation given a rotation matrix 
#### Tsagris Michail 11/2013 
#### mtsagris@yahoo.gr
#### References: Howard E Haber couse webpage
#### http://scipp.ucsc.edu/~haber/ph216/rotation_12.pdf 
################################

Arotation <- function(A) {
  ## A is a 3x3 rotation matrix
  con1 = round(det(A), 15)
  con2 = round( mean( abs( A %*% t(A) - diag(3) ) ), 15 )
  if ( con1 != 1 |  con2 > .Machine$double.eps ) {
    res <- paste("This is not a rotation matrix")
  } else {
    tr <- sum(diag(A))
    rad <- acos(0.5 * (tr - 1))
    angle <- rad * 180 / pi  ## from rads to degrees
    ksi <- c(A[3, 2] - A[2, 3], A[1, 3] - A[3, 1], A[2, 1] - A[1, 2])/
    sqrt((3 - tr) * (1 + tr))
    axis <- c( asin(ksi[3]), atan2(ksi[2], ksi[1]) )
    axis <- c(axis / pi * 180)  ## from degrees to rads
    ## if the latitude or longitude are negative add 360 (degrees) 
    axis[axis<0] <- axis[axis<0] + 360
    names(axis) <- c("latitude", "longitude")
    res <- list(angle = angle, axis = axis)
  }
  res
}
