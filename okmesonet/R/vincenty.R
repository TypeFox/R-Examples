vincenty <- function(lon1, lat1, lon2, lat2) {
  ## Function to return the Vincenty distance for two points, using WGS84 
  ## ellipsoid
  ##
  ## code adapted from JavaScript written by Chris Veness
  ## http://www.movable-type.co.uk/scripts/latlong-vincenty.html
  ##
  ## Vincenty Inverse Solution of Geodesics on the Ellipsoid 
  ## (c) Chris Veness 2002-2011 
  ## from: Vincenty inverse formula - T Vincenty, "Direct and Inverse 
  ##  Solutions of Geodesics on the Ellipsoid with application of nested 
  ##  equations", Survey Review, vol XXII no 176, 1975
  ##  http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
  ##
  ## Arguments:
  ##  lon1: longitude of point 1, decimal degrees
  ##  lat1: latitude of point 1, decimal degrees
  ##  lon2: longitude of point 2, decimal degrees
  ##  lat1: latitude of point 1, decimal degrees
  ## Returns: great circle distance in meters between two points
  
  ## WGS84 ellipsoid parameters 
  a <- 6378137
  b <- 6356752.314245
  f <- 1/298.257223563
  
  ## convert to radians
  lon1 <- lon1 * (pi/180)
  lat1 <- lat1 * (pi/180)
  lon2 <- lon2 * (pi/180)
  lat2 <- lat2 * (pi/180)
  
  L <- (lon2 - lon1)
  U1 <- atan((1 - f) * tan(lat1))
  U2 <- atan((1 - f) * tan(lat2))
  sinU1 <- sin(U1)
  cosU1 <- cos(U1)
  sinU2 <- sin(U2)
  cosU2 <- cos(U2)
  lambda <- L
  lambdaP <- runif(1, 1, 2) # random number to begin while loop
  
  iterLimit <- 0
  while ((abs(lambda-lambdaP) > 1e-12) && iterLimit < 100) {
    sinLambda <- sin(lambda)
    cosLambda <- cos(lambda)
    sinSigma <- sqrt((cosU2 * sinLambda) * (cosU2 * sinLambda) 
                     + (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) 
                     * (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda))
    if (sinSigma==0) break # identical points
    cosSigma <- sinU1 * sinU2 + cosU1 * cosU2 * cosLambda
    sigma <- atan2(sinSigma, cosSigma)
    sinAlpha <- cosU1 * cosU2 * sinLambda / sinSigma
    cosSqAlpha <- 1 - sinAlpha * sinAlpha
    cos2SigmaM <- cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha    
    if (is.nan(cos2SigmaM)) cos2SigmaM <- 0 
    C <- f/16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha))
    lambdaP <- lambda
    lambda <- L + (1 - C) * f * sinAlpha * (sigma + C * sinSigma * 
      (cos2SigmaM + C * cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM)))
    iterLimit <- iterLimit + 1
  }
  
  ## return 0 if identical points
  if (sinSigma==0) {
    return(0)
  } else if (iterLimit==100) {
    return(NA)  # failed to converge; return NA
  } else {
    uSq <- cosSqAlpha * (a * a - b * b) / (b * b)
    A <- 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
    B <- uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))
    deltaSigma <- B * sinSigma * 
      (cos2SigmaM + B / 4 * (cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM) 
                             - B / 6 * cos2SigmaM 
                             * (-3 + 4 * sinSigma * sinSigma) 
                             * (-3 + 4 * cos2SigmaM * cos2SigmaM)))
    s  <- b * A * (sigma - deltaSigma)
    return(round(s, digits=0)) # round to meter precision
  }
}