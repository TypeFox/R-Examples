# ------------------------------------------------------------------------------
# Estimating the surface of 1 dg by 1 dg cells on earth
# ------------------------------------------------------------------------------
# based on data that give the surface distance per 1 dg change in lat/lon
# Data below are from http://en.wikipedia.org/wiki/Latitude
# Latdist, Londist =  Surface distance per 1 deg change in latitude/longitude
# Latrad = N-S radius of curvature, Lat, Lonrad = E-W radius


earth_surf <- function (lat = 0, lon = 0) {
  Earth <- data.frame(lat = seq(0, 90, by = 15),
    Latdist = c(110.574, 110.649, 110.852, 111.132, 111.412, 111.618, 111.694),
    Londist = c( 111.320, 107.551, 96.486, 78.847, 55.800 , 28.902, 0)        ,
    Latrad  = c(6335.44, 6339.70, 6351.38, 6367.38, 6383.45, 6395.26, 6399.59) ,
    Lonrad  = c(6378.14, 6379.57, 6383.48, 6388.84, 6394.21, 6398.15, 6399.59) )

  Latdist <- approx(Earth$lat, Earth$Latdist, abs(lat))$y
  Londist <- approx(Earth$lat, Earth$Londist, abs(lat))$y
  Latdist*Londist * 1.008539 *1e6 #correct ~ 1 % error !
}

earth_dist <- function (alat, alon, blat, blon, method = 1)   {

  if (any(alat > 90)  | any(blat > 90) | any(alat <  -90) | any(blat <  -90))
    stop("'alat' and 'blat' should be inbetween -90 and 90")
  if (any(alon > 180) | any(blon > 180)| any(alon < -180) | any(blon < -180))
    stop("'alat' and 'blat' should be inbetween -90 and 90")

# Distance in meters between two latitude/longitude pairs
# assumes the earth is a perfect sphere
  radius <- 6371000 # earth radius in metres

  alon <- alon * pi / 180
  alat <- alat * pi / 180
  blon <- blon * pi / 180
  blat <- blat * pi / 180

  if (method == 1) {
    fac <- cos(alat) * cos(blat) * (cos(alon) * cos(blon) + sin(alon)* sin(blon)) +
       sin(alat) * sin(blat)
    DD <- acos(fac) * radius
  } else {
    dLat = (blat-alat)
    dLon = (blon-alon)
    A   <- sin(dLat/2) * dLat/2 + sin(dLon/2) * dLon/2 * cos(alat) * cos(blat)
    fac <- 2 * atan2(sqrt(A), sqrt(1.-A)) # angular distance in radians
    DD  <- radius * fac
  }
  return(DD)
}


