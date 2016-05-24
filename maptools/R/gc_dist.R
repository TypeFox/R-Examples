# 2007 Eric Archer and Roger Bivand
# 

gcDestination <- function(lon, lat, bearing, dist, dist.units = "km",
    model=NULL, Vincenty=FALSE) {
 # lat, lon : lattitude and longitude in decimal degrees
 # bearing : bearing from 0 to 360 degrees
 # dist : distance travelled
 # dist.units : units of distance "km" (kilometers), "nm" (nautical 
 # miles), "mi" (statute miles)
 # model : choice of ellipsoid model ("WGS84", "GRS80", "Airy", 
 # "International", "Clarke", "GRS67")

    if (!is.numeric(lon)) stop("lon not numeric")
    if (!is.numeric(lat)) stop("lat not numeric")
    if (!is.numeric(bearing)) stop("bearing not numeric")
    if (!is.numeric(dist)) stop("dist not numeric")

    if (length(lon) != length(lat)) stop("lon and lat differ in length")
    if (length(bearing) > 1L && length(lon) > 1L) stop("length mismatch")
    if (length(bearing) > 1L && length(dist) > 1L) stop("length mismatch")

    as.radians <- function(degrees) degrees * pi / 180
    as.degrees <- function(radians) radians * 180 / pi
    as.bearing <- function(radians) (as.degrees(radians) + 360) %% 360

    ellipsoid <- function(model = "WGS84") {
        switch(model,
        WGS84 = c(a = 6378137, b = 6356752.3142, f = 1 / 298.257223563),
        GRS80 = c(a = 6378137, b = 6356752.3141, f = 1 / 298.257222101),
        Airy = c(a = 6377563.396, b = 6356256.909, f = 1 / 299.3249646),
        International = c(a = 6378888, b = 6356911.946, f = 1 / 297),
        Clarke = c(a = 6378249.145, b = 6356514.86955, f = 1 / 293.465),
        GRS67 = c(a = 6378160, b = 6356774.719, f = 1 / 298.25),
        c(a = NA, b = NA, f = NA)
    )}
    
    dist <- switch(dist.units,
        km = dist,
        nm = dist * 1.852,
        mi = dist * 1.609344
    )
    lat <- as.radians(lat)
    lon <- as.radians(lon)
    bearing <- as.radians(bearing)

    if (is.null(model)) {
 # Code adapted from JavaScript by Chris Veness 
 # (scripts@movable-type.co.uk) at
 # http://www.movable-type.co.uk/scripts/latlong.html#ellipsoid
 #   originally from Ed Williams' Aviation Formulary, 
 # http://williams.best.vwh.net/avform.htm
        radius <- 6371
        psi <- dist / radius
        lat2 <- asin(sin(lat) * cos(psi) +  cos(lat) * sin(psi) * cos(bearing))
        lon2 <- lon + atan2(sin(bearing) * sin(psi) * cos(lat), cos(psi) - 
            sin(lat) * sin(lat2))
        if (is.nan(lat2) || is.nan(lon2)) warning("Out of range values")
        return(cbind(long=as.degrees(lon2), lat=as.degrees(lat2)))
    }

    ellips <- ellipsoid(model)
    if (is.na(ellips["a"])) stop("no such ellipsoid model")
    if (Vincenty) {
 # Code adapted from JavaScript by Chris Veness 
 # (scripts@movable-type.co.uk) at
 # http://www.movable-type.co.uk/scripts/latlong-vincenty-direct.html
 # Original reference (http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf):
 #   Vincenty, T. 1975.  Direct and inverse solutions of geodesics on 
 # the ellipsoid with application of nested equations.
 #      Survey Review 22(176):88-93
        dist <- dist * 1000
        sin.alpha1 <- sin(bearing)
        cos.alpha1 <- cos(bearing)
        tan.u1 <- (1 - ellips["f"]) * tan(lat)
        cos.u1 <- 1 / sqrt(1 + (tan.u1 ^ 2))
        sin.u1 <- tan.u1 * cos.u1
        sigma1 <- atan2(tan.u1, cos.alpha1)
        sin.alpha <- cos.u1 * sin.alpha1
        cos.sq.alpha <- 1 - (sin.alpha ^ 2)
        u.sq <- cos.sq.alpha * ((ellips["a"] ^ 2) - (ellips["b"] ^ 2)) /
            (ellips["b"] ^ 2)
        cap.A <- 1 + u.sq / 16384 * (4096 + u.sq * (-768 + u.sq * (320 - 
            175 * u.sq)))
        cap.B <- u.sq / 1024 * (256 + u.sq * (-128 + u.sq * (74 - 47 * u.sq)))

        sigma <- dist / (ellips["b"] * cap.A)
        sigma.p <- 2 * pi
        cos.2.sigma.m <- cos(2 * sigma1 + sigma)
        while(any(abs(sigma - sigma.p) > 1e-12)) {
            cos.2.sigma.m <- cos(2 * sigma1 + sigma)
            sin.sigma <- sin(sigma)
            cos.sigma <- cos(sigma)
            delta.sigma <- cap.B * sin.sigma * (cos.2.sigma.m + cap.B / 4 * 
                (cos.sigma *
                (-1 + 2 * cos.2.sigma.m ^ 2) - cap.B / 6 * cos.2.sigma.m *
                (-3 + 4 * sin.sigma ^ 2) * (-3 + 4 * cos.2.sigma.m ^ 2)))
            sigma.p <- sigma
            sigma <- dist / (ellips["a"] * cap.A) + delta.sigma
        }
        tmp <- sin.u1 * sin.sigma - cos.u1 * cos.sigma * cos.alpha1
        lat2 <- atan2(sin.u1 * cos.sigma + cos.u1 * sin.sigma * cos.alpha1,
            (1 - ellips["f"]) * sqrt(sin.alpha ^ 2 + tmp ^ 2))
        lambda <- atan2(sin.sigma * sin.alpha1, cos.u1 * cos.sigma - sin.u1 * 
            sin.sigma * cos.alpha1)
        cap.C <- ellips["f"] / 16 * cos.sq.alpha * (4 + ellips["f"] * 
            (ellips["f"] - 3 * cos.sq.alpha))
        cap.L <- lambda - (1 - cap.C) * ellips["f"] * sin.alpha *
            (sigma + cap.C * sin.sigma * (cos.2.sigma.m + cap.C * cos.sigma * 
            (-1 + 2 * cos.2.sigma.m ^ 2)))
        lat2 <- as.degrees(lat2)
        lon2 <- as.degrees(lon + cap.L)
    } else {
 # Code adapted from JavaScript by Larry Bogan (larry@go.ednet.ns.ca) 
 # at http://www.go.ednet.ns.ca/~larry/bsc/jslatlng.html
        e <- 0.08181922
        radius <- (ellips["a"] / 1000) * (1 - e^2) / ((1 - e^2 * 
            sin(lat)^2)^1.5)
        psi <- dist / radius
        phi <- pi / 2 - lat
        arc.cos <- cos(psi) * cos(phi) + sin(psi) * sin(phi) * cos(bearing)
        lat2 <- as.degrees((pi / 2) - acos(arc.cos))
        arc.sin <- sin(bearing) * sin(psi) / sin(phi)
        lon2 <- as.degrees(lon + asin(arc.sin))
    }
    return(cbind(long=lon2, lat=lat2))
}
