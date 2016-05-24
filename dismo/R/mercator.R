# Author: Robert J. Hijmans
# Date :  July 2010
# Version 1.0
# Licence GPL v3

# for use with gmap

Mercator <- function (p, inverse = FALSE) {
	r = 6378137
    toRad <- pi/180
    if (inverse) {
        p <- .pointsToMatrix(p, checkLonLat = FALSE)
        p[, 2] <- pi/2 - 2 * atan(exp(-p[, 2]/r))
        p[, 1] <- p[, 1]/r
        colnames(p) <- c("lon", "lat")
        return(p/toRad)
    }
    else {
        p <- .pointsToMatrix(p) * toRad
        p[, 2] <- log(tan(p[, 2]) + (1/cos(p[, 2])))
        p <- p * r
        colnames(p) <- c("x", "y")
        return(p)
    }
}

