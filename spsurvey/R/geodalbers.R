geodalbers <- function(lon, lat, sph="GRS80", clon=-96, clat=23, sp1=29.5,
   sp2=45.5) {

################################################################################
# Function: geodalbers
# Purpose: Project latitude and longitude (spheroid) to Albers projection in
#   plane
# Programmer: Tony Olsen translated from Denis White C program
# Date: September 4, 2003 (White November 2001)
# Arguments:
#	  lon = longitude (decimal degrees) vector to be projected using Albers.
# 	lat = latitude (decimal degrees) vector to be projected using Albers.
#   sph = Spheroid options: Clarke1866, GRS80, WGS84.  The default is
#     GRS80.
#   clon = Center longitude (decimal degrees).  The default is -96.
#   clat = Origin latitude (decimal degrees).  The default is 23.
#   sp1 = Standard parallel 1 (decimal degrees).  The default is 29.5.
#   sp2 = Standard parallel 2 (decimal degrees).  The default is 45.5.
# Output:
#   A data frame of Albers x-coordinate and y-coordinate projections for
#   latitude and longitude.
################################################################################

# Specify parameters for selected spheroid
  if(sph == "Clarke1866") { 
    a <- 6378206.4
    b <- 6356583.8 
  } else if(sph == "GRS80") { 
    a <- 6378137.0
    b <- 6356752.31414
  } else if(sph == "WGS84") {
    a <- 6378137.0
    b <- 6356752.31424518
  } else {
    stop("\nSpheroid does not match available options.")
  }

# Do calculations for selected spheroid
  DEGRAD <- (pi/180)
  clat <- clat*DEGRAD 
  clon <- clon*DEGRAD
  sp1 <- sp1*DEGRAD
  sp2 <- sp2*DEGRAD
  e2 <- 1.0 - (b*b) / (a*a)
  e4 <- e2 * e2 
  e <- sqrt(e2)
  t1 <- 1.0 - e2
  t2 <- 1.0 / (2.0*e)
  sinlat <- sin(clat) 
  t3 <- 1.0 - e2*sinlat*sinlat
  q0 <- t2 * log((1.0 - e*sinlat)/(1.0 + e*sinlat))
  q0 <- t1 * (sinlat/t3 - q0)
  sinlat <- sin(sp1) 
  t3 <- 1.0 - e2*sinlat*sinlat
  q1 <- t2 * log((1.0 - e*sinlat)/(1.0 + e*sinlat))
  q1 <- t1 * (sinlat/t3 - q1)
  m1 <- cos(sp1) / sqrt(t3)
  sinlat <- sin(sp2) 
  t3 <- 1.0 - e2*sinlat*sinlat
  q2 <- t2 * log((1.0 - e*sinlat)/(1.0 + e*sinlat))
  q2 <- t1 * (sinlat/t3 - q2)
  m2 <- cos(sp2) / sqrt(t3)
  n <- (m1*m1 - m2*m2) / (q2 - q1)
  C <- m1*m1 + n*q1
  rho0 <- a * sqrt(C - n*q0) / n

# project from lon/lat to plane with Albers conic projection
  lat <- lat*DEGRAD
  lon <- lon*DEGRAD
  sinlat <- sin(lat) 
  q <- sinlat / (1.0 - e2 * sinlat*sinlat)
  q <- t1 * (q - t2 * log ((1.0 - e*sinlat)/(1.0 + e*sinlat)))
  rho <- a * sqrt(C - n*q) / n
  theta <- n * (lon - clon)
  xyalbers <- data.frame(xcoord=rho*sin(theta), ycoord=rho0 - rho*cos(theta))

# Return the result
   xyalbers
}
