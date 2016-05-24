# Copyright (c) 2006 Roger Bivand
qibla <- function(from, type="abdali") {
# Mecca as given by K. Abdali
  to <- c(39.82333, 21.42333)
  gzAzimuth(from=from, to=to, type=type)
}
# http://www.geocities.com/sualeh85/sfweb/QTable.html
# http://patriot.net/users/abdali/ftp/qibla.pdf
# http://www.patriot.net/users/abdali/ftp/praytimer.zip
# http://www.world-gazetteer.com/

gzAzimuth <- function(from, to, type="snyder_sphere") {
  deg2rad <- function(x) x*pi/180
  rad2deg <- function(x) x*180/pi
# note negative longitudes
  if (is.matrix(from)) {
    lon <- -deg2rad(from[,1])
    lat <- deg2rad(from[,2])
  } else {
    lon <- -deg2rad(from[1])
    lat <- deg2rad(from[2])
  }
  if (is.matrix(to)) {
    if (nrow(to) > 1) stop("to: single coordinate")
    to <- c(to)
  } 
  lon0 <- -deg2rad(to[1])
  lat0 <- deg2rad(to[2])
# bug found by Sebastian Luque
  dflon = lon-lon0
# results in degrees from N, negative west
  if (type == "abdali") res <- atan2(sin(dflon), ((cos(lat)*tan(lat0)) - 
    (sin(lat)*cos(dflon))))
  else if (type == "snyder_sphere") res <- atan2((cos(lat0)*sin(dflon)), 
    (cos(lat)*sin(lat0)) - (sin(lat)*cos(lat0)*cos(dflon)))
  else stop("type unkown")
  is.na(res) <- lon == lon0 & lat == lat0 
  rad2deg(res)
}

trackAzimuth <- function(track, type="snyder_sphere") {
  if (!is.matrix(track)) stop("track must be two-column matrix")
  if (ncol(track) != 2) stop("track must be two-column matrix")
  n1 <- nrow(track)-1
  if (n1 < 2) stop("less than two points")
  res <- numeric(n1)
  for (i in seq(along=res)) res[i] <- gzAzimuth(track[i,], track[(i+1),],
    type=type)
  res
} 


