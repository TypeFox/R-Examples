#
### Great circle distances among sites
### Original code by Mario Pineda-Krch, taken from
### http://www.r-bloggers.com/great-circle-distance-calculations-in-r/
### and adaptated.
#
# Geodesic distance using the Spherical Law of Cosines (slc)
gcd.slc <- function(ss,radius=6371) {
  rr <- matrix(NA,nrow(ss),nrow(ss),dimnames=list(rownames(ss),rownames(ss)))
  for(j in 1L:(nrow(ss)-1L)) {
    rr[(j+1L):nrow(ss),j] <- acos(sin(ss[j,1L]*pi/180)*sin(ss[(j+1L):nrow(ss),1L]*pi/180)+
                                  cos(ss[j,1L]*pi/180)*cos(ss[(j+1L):nrow(ss),1L]*pi/180)*
                                  cos((ss[(j+1L):nrow(ss),2L]-ss[j,2L])*pi/180)) * radius
  }
  return(as.dist(rr))
}
#
# Geodesic distance using the Haversine formula (hf)
gcd.hf <- function(ss,radius=6371) {
  rr <- matrix(NA,nrow(ss),nrow(ss),dimnames=list(rownames(ss),rownames(ss)))
  for(j in 1L:(nrow(ss)-1L)) {
    delta.long <- (ss[(j+1L):nrow(ss),2L]-ss[j,2L])*pi/180
    delta.lat <- (ss[(j+1L):nrow(ss),1L]-ss[j,1L])*pi/180
    a <- sin(delta.lat/2)^2 + cos(ss[j,1L]*pi/180) * cos(ss[(j+1L):nrow(ss),1L]*pi/180) * sin(delta.long/2)^2
    d <- numeric(length(a))
    for(i in 1L:length(a)) d[i] <- 2 * asin(min(1,sqrt(a[i]))) * radius
    rr[(j+1L):nrow(ss),j] <- d
  }
  return(as.dist(rr))
}
#
# Calculate single geodesic distance between two points using Vincenty inverse formula for ellipsoids (vife)
.vife <- function(long1, lat1, long2, lat2, a, b, f, iterLimit = 100) {
  L <- long2-long1
  U1 <- atan((1-f) * tan(lat1))
  U2 <- atan((1-f) * tan(lat2))
  sinU1 <- sin(U1)
  cosU1 <- cos(U1)
  sinU2 <- sin(U2)
  cosU2 <- cos(U2)
  cosSqAlpha <- NULL
  sinSigma <- NULL
  cosSigma <- NULL
  cos2SigmaM <- NULL
  sigma <- NULL
  lambda <- L
  lambdaP <- 0
  while (abs(lambda-lambdaP) > 1e-12 & iterLimit>0) {
    sinLambda <- sin(lambda)
    cosLambda <- cos(lambda)
    sinSigma <- sqrt( (cosU2*sinLambda) * (cosU2*sinLambda) +
                      (cosU1*sinU2-sinU1*cosU2*cosLambda) * (cosU1*sinU2-sinU1*cosU2*cosLambda) )
    if (sinSigma==0) return(0)  # Co-incident points
    cosSigma <- sinU1*sinU2 + cosU1*cosU2*cosLambda
    sigma <- atan2(sinSigma, cosSigma)
    sinAlpha <- cosU1 * cosU2 * sinLambda / sinSigma
    cosSqAlpha <- 1 - sinAlpha*sinAlpha
    cos2SigmaM <- cosSigma - 2*sinU1*sinU2/cosSqAlpha
    if (is.na(cos2SigmaM)) cos2SigmaM <- 0  # Equatorial line: cosSqAlpha=0
    C <- f/16*cosSqAlpha*(4+f*(4-3*cosSqAlpha))
    lambdaP <- lambda
    lambda <- L + (1-C) * f * sinAlpha *
              (sigma + C*sinSigma*(cos2SigmaM+C*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)))
    iterLimit <- iterLimit - 1
  }
  if (iterLimit==0) stop("Failed to converge")
  uSq <- cosSqAlpha * (a*a - b*b) / (b*b)
  A <- 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)))
  B <- uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)))
  deltaSigma = B*sinSigma*(cos2SigmaM+B/4*(cosSigma*(-1+2*cos2SigmaM^2) -
                                      B/6*cos2SigmaM*(-3+4*sinSigma^2)*(-3+4*cos2SigmaM^2)))
  s <- b*A*(sigma-deltaSigma) / 1000
  return(s) # Distance in km
}
#
# Geodesic distance using the Vincenty inverse formula for ellipsoids (vife)
gcd.vife <- function(ss, a = 6378137, b = 6356752.314245, f = 1/298.257223563) {
  rr <- matrix(NA,nrow(ss),nrow(ss),dimnames=list(rownames(ss),rownames(ss)))
  for(j in 1L:(nrow(ss)-1L)) {
    for(i in (j+1L):nrow(ss)) {
      if(ss[j,2L] == ss[i,2L]) {
        rr[i,j] <- NA
        warning(paste("Location ",i," have the same longitude as location ",j,
                      ", the Vincenty inverse formula can therefore not be applied. ",
                      "Use an alternative approach.",sep=""))
      } else {
        rr[i,j] <- .vife(ss[j,2L]*pi/180, ss[j,1L]*pi/180, ss[i,2L]*pi/180, ss[i,1L]*pi/180, a, b, f)
      }
    }
  }
  return(as.dist(rr))
}
#
