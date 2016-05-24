##' @title Sinusoidal projection
##' @param r Lattitude-longitude coordinates in a matrix with columns
##' labelled \code{phi} (lattitude) and \code{lambda}
##' (longitude). Alternatively string "boundary", indicating that
##' boundary of projection should be drawn.
##' @param proj.centre Location of centre of projection as matrix with
##' column names \code{phi} (elevation) and \code{lambda}
##' (longitude). Currently only longitude is used by this function.
##' @param lambdalim Limits of longitude to plot
##' @param lines If this is \code{TRUE} create breaks of \code{NA}s
##' when lines cross the limits of longitude. This prevents lines
##' crossing the centre of the projection.
##' @param ... Arguments not used by this projection.
##' @return Two-column matrix with columns labelled \code{x} and
##' \code{y} of locations of projection of coordinates on plane
##' @references \url{http://en.wikipedia.org/wiki/Map_projection},
##' \url{http://mathworld.wolfram.com/SinusoidalProjection.html}
##' @author David Sterratt
##' @export
sinusoidal <- function(r, proj.centre=cbind(phi=0, lambda=0),
                       lambdalim=NULL, lines=FALSE, ...) {
  if (is.character(r) & identical(r, "boundary")) {
    return(cbind(x=NA, y=NA))
  }
  lambda0 <- proj.centre[1, "lambda"]
  x <- (r[,"lambda"] - lambda0)*cos(r[,"phi"])
  y <- r[,"phi"]
  rc <- cbind(x=x, y=y)
  if (!is.null(lambdalim)) {
    rc[r[,"lambda"] > lambdalim[2],] <- NA
    rc[r[,"lambda"] < lambdalim[1],] <- NA
    if (lines) {
      inds <- which(abs(diff(r[,"lambda"])) > 0.5*(lambdalim[2] - lambdalim[1]))
      rc[inds,] <- NA
    }
  }
  return(rc)
}

##' @title Orthographic projection
##' @param r Lattitude-longitude coordinates in a matrix with columns
##' labelled \code{phi} (lattitude) and \code{lambda} (longitude)
##' @param proj.centre Location of centre of projection as matrix with
##' column names \code{phi} (elevation) and \code{lambda} (longitude).
##' @param ... Arguments not used by this projectio.n
##' @return Two-column matrix with columns labelled \code{x} and
##' \code{y} of locations of projection of coordinates on plane
##' @references \url{http://en.wikipedia.org/wiki/Map_projection},
##' \url{http://mathworld.wolfram.com/OrthographicProjection.html}
##' @author David Sterratt
##' @export
orthographic <- function(r,
                         proj.centre=cbind(phi=0, lambda=0),
                         ...) {
  if (is.character(r) & identical(r, "boundary")) {
    return(circle(360))
  }
  lambda0 <- proj.centre[1, "lambda"]
  phi0    <- proj.centre[1, "phi"]
  
  ## First translate to Cartesian coordinates with only the rotation
  ## about the polar axis so that (0, lambda0) is at the centre of the
  ## projection.
  P <- cbind(x=cos(r[,"phi"])*sin(r[,"lambda"] - lambda0),
             y=sin(r[,"phi"]),
             z=cos(r[,"phi"])*cos(r[,"lambda"] - lambda0))

  ## The rotate about the x axis so that (phi0, lambda0) is at the
  ## centre
  P <- P %*% rbind(c(1, 0, 0),
                   c(0,  cos(phi0), sin(phi0)),
                   c(0, -sin(phi0), cos(phi0)))
  
  ## Projection onto the x-y plane
  rc <- cbind(x=P[, 1], y=P[, 2])

  ## Anything with negative z is obsured
  rc[P[, 3] < 0,] <- NA
  return(rc)
}

##' @title Lambert azimuthal equal area projection
##' @param r 2-column Matrix of spherical coordinates of points on
##' sphere. Column names are \code{phi} and \code{lambda}.
##' @param ... Arguments not used by this projection.
##' @return 2-column Matrix of Cartesian coordinates of points on polar
##' projection. Column names should be \code{x} and \code{y}.
##' @author David Sterratt
##' @note This is a special case with the point centred on the
##' projection being the South Pole. The Mathworld equations are for
##' the more general case.
##' @references \url{http://en.wikipedia.org/wiki/Map_projection},
##' \url{http://mathworld.wolfram.com/LambertAzimuthalEqual-AreaProjection.html}
##' Fisher, N. I., Lewis, T., and Embleton,
##' B. J. J. (1987). Statistical analysis of spherical data. Cambridge
##' University Press, Cambridge, UK.
##' @export
azimuthal.equalarea <- function(r, ...) {
  if (is.character(r) & identical(r, "boundary")) {
    return(cbind(x=NA, y=NA))
  }
  rho <- sqrt(2*(1 + sin(r[,"phi"])))
  x <- rho*cos(r[,"lambda"])
  y <- rho*sin(r[,"lambda"])
  return(cbind(x=x, y=y))
}

##' @title Azimuthal equidistant projection
##' @param r 2-column Matrix of spherical coordinates of points on
##' sphere. Column names are \code{phi} and \code{lambda}.
##' @param ... Arguments not used by this projection.
##' @return 2-column Matrix of Cartesian coordinates of points on polar
##' projection. Column names should be \code{x} and \code{y}.
##' @author David Sterratt
##' @note This is a special case with the point centred on the
##' projection being the South Pole. The Mathworld equations are for
##' the more general case.
##' @references \url{http://en.wikipedia.org/wiki/Map_projection},
##' \url{http://mathworld.wolfram.com/AzimuthalEquidistantProjection.html}  
##' @export
azimuthal.equidistant <- function(r, ...) {
  if (is.character(r) & identical(r, "boundary")) {
    return(cbind(x=NA, y=NA))
  }
  rho <- pi/2 + r[,"phi"]
  x <- rho*cos(r[,"lambda"])
  y <- rho*sin(r[,"lambda"])
  return(cbind(x=x, y=y))
}

##' @title Azimuthal conformal or stereographic or Wulff projection
##' @param r 2-column Matrix of spherical coordinates of points on
##' sphere. Column names are \code{phi} and \code{lambda}.
##' @param ... Arguments not used by this projection.
##' @return 2-column Matrix of Cartesian coordinates of points on polar
##' projection. Column names should be \code{x} and \code{y}.
##' @author David Sterratt
##' @note This is a special case with the point centred on the
##' projection being the South Pole. The Mathworld equations are for
##' the more general case.
##' @references \url{http://en.wikipedia.org/wiki/Map_projection},
##' \url{http://mathworld.wolfram.com/StereographicProjection.html}
##' Fisher, N. I., Lewis, T., and Embleton,
##' B. J. J. (1987). Statistical analysis of spherical data. Cambridge
##' University Press, Cambridge, UK.
##' @export
azimuthal.conformal <- function(r, ...) {
  if (is.character(r) & identical(r, "boundary")) {
    return(cbind(x=NA, y=NA))
  }
  rho <- tan(pi/4 + r[,"phi"]/2)
  x <- rho*cos(r[,"lambda"])
  y <- rho*sin(r[,"lambda"])
  return(cbind(x=x, y=y))
}

