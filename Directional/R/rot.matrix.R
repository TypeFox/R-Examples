################################
#### Rotation matrix from a rotation axis and angle of rotation
#### Tsagris Michail 11/2013 
#### mtsagris@yahoo.gr
#### References: Chang Ted (1986)
#### Spherical egession. Annals of statistics, 14(3): 907-924
################################

rot.matrix <- function(ksi, theta, rads = FALSE) {
  ## ksi is the rotation axis, where the first element is the
  ## latitude and the second is the longitude
  ## theta is the angle of rotation
  if (rads == TRUE) {
    lat <- ksi[1]
    long <- ksi[2]
    the <- theta
  } else {
    lat <- ksi[1] * pi/180
    long <- ksi[2] * pi/180
    the <- theta * pi/180
  }
  t1 <- cos(lat) * cos(long)
  t2 <- cos(lat) * sin(long)
  t3 <- sin(lat)
  L <- matrix(c(0, t3, -t2, -t3, 0, t1, t2, -t1, 0), ncol = 3)
  diag(3) + L * sin(the) + crossprod(L) * ( 1 - cos(the) )
}