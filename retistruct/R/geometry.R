##
## Geometry functions
## 

##' @title Vector norm
##' @param X Vector or matrix. 
##' @return If a vector, returns the 2-norm  of the
##' vector. If a matrix, returns the 2-norm of each row of the matrix
##' @author David Sterratt
##' @export
vecnorm <- function(X) {
  if (is.vector(X)) {
    return(sqrt(sum(X^2)))
  } else {
    return(sqrt(rowSums(X^2)))
  }
}

##' @title "Signed area" of triangles on a plane
##' @param P 2-column matrix of vertices of triangles
##' @param Pt 3-column matrix of indicies of rows of \code{P} giving
##' triangulation
##' @return Vectors of signed areas of triangles. Positive sign
##' indicates points are anticlockwise direction; negative indicates
##' clockwise.
##' @author David Sterratt
##' @export
tri.area.signed <- function(P, Pt) {
  A <- P[Pt[,1],]
  B <- P[Pt[,2],]
  C <- P[Pt[,3],]
  AB <- cbind(B-A, 0)
  BC <- cbind(C-B, 0)
  return(0.5 * extprod3d(AB, BC)[,3])
}

##' @title Area of triangles on a plane
##' @param P 2-column matrix of vertices of triangles
##' @param Pt 3-column matrix of indicies of rows of \code{P} giving
##' triangulation
##' @return Vectors of areas of triangles
##' @author David Sterratt
##' @export
tri.area <- function(P, Pt) {
  return(abs(tri.area.signed(P, Pt)))
}

##' This uses l'Hullier's theorem to compute the spherical excess and
##' hence the area of the spherical triangle.
##' 
##' @title Area of triangles on a sphere
##' @param P 2-column matrix of vertices of triangles given in
##' spherical polar coordinates. Columns need to be labelled
##' \code{phi} (lattidute) and \code{lambda} (longitude).
##' @param Pt 3-column matrix of indicies of rows of \code{P} giving
##' triangulation
##' @return Vectors of areas of triangles in units of steradians
##' @source Wolfram MathWorld
##' \url{http://mathworld.wolfram.com/SphericalTriangle.html} and
##' \url{http://mathworld.wolfram.com/SphericalExcess.html}
##' @author David Sterratt
##' @examples
##' ## Something that should be an eighth of a sphere, i.e. pi/2
##' P <- cbind(phi=c(0, 0, pi/2), lambda=c(0, pi/2, pi/2))
##' Pt <- cbind(1, 2, 3)
##' ## The result of this should be 0.5
##' print(sphere.tri.area(P, Pt)/pi)
##'
##' ## Now a small triangle
##' P1 <- cbind(phi=c(0, 0, 0.01), lambda=c(0, 0.01, 0.01))
##' Pt1 <- cbind(1, 2, 3)
##' ## The result of this should approximately 0.01^2/2
##' print(sphere.tri.area(P, Pt)/(0.01^2/2))
##'
##' ## Now check that it works for both 
##' P <- rbind(P, P1)
##' Pt <- rbind(1:3, 4:6)
##' ## Should have two components
##' print(sphere.tri.area(P, Pt))
##' @export
sphere.tri.area <- function(P, Pt) {
  ## Get lengths of sides of all triangles
  a <- central.angle(P[Pt[,1],"phi"], P[Pt[,1],"lambda"],
                     P[Pt[,2],"phi"], P[Pt[,2],"lambda"])
  b <- central.angle(P[Pt[,2],"phi"], P[Pt[,2],"lambda"],
                     P[Pt[,3],"phi"], P[Pt[,3],"lambda"])
  c <- central.angle(P[Pt[,3],"phi"], P[Pt[,3],"lambda"],
                     P[Pt[,1],"phi"], P[Pt[,1],"lambda"])
  
  ## Semiperimeter is half perimeter
  s <- 1/2*(a + b + c)

  ## Compute excess using L'Huilier's Theorem
  E <- 4*atan(sqrt(tan(s/2)*tan((s-a)/2)*tan((s-b)/2)*tan((s-c)/2)))
  names(E) <- c()
  return(E)
}

##' Determine the intersection of two lines L1 and L2 in two dimensions,
##' using the formula described by Weisstein.
##' 
##' @title Determine intersection between two lines 
##' @param P1 vector containing x,y coordinates of one end of L1
##' @param P2 vector containing x,y coordinates of other end of L1
##' @param P3 vector containing x,y coordinates of one end of L2
##' @param P4 vector containing x,y coordinates of other end of L2
##' @param interior.only boolean flag indicating whether only
##' intersections inside L1 and L2 should be returned.
##' @return Vector containing x,y coordinates of intersection of L1
##' and L2.  If L1 and L2 are parallel, this is infinite-valued.  If
##' \code{interior.only} is \code{TRUE}, then when the intersection
##' does not occur between P1 and P2 and P3 and P4, a vector
##' containing \code{NA}s is returned.
##' @source Weisstein, Eric W. "Line-Line Intersection."
##' From MathWorld--A Wolfram Web Resource.
##' \url{http://mathworld.wolfram.com/Line-LineIntersection.html}
##' @author David Sterratt
##' @export
##' @examples
##' ## Intersection of two intersecting lines
##' line.line.intersection(c(0, 0), c(1, 1), c(0, 1), c(1, 0))
##'
##' ## Two lines that don't intersect
##' line.line.intersection(c(0, 0), c(0, 1), c(1, 0), c(1, 1))
line.line.intersection <- function(P1, P2, P3, P4, interior.only=FALSE) {
  P1 <- as.vector(P1)
  P2 <- as.vector(P2)
  P3 <- as.vector(P3)
  P4 <- as.vector(P4)

  dx1 <- P1[1] - P2[1]
  dx2 <- P3[1] - P4[1]
  dy1 <- P1[2] - P2[2]
  dy2 <- P3[2] - P4[2]

  D <- det(rbind(c(dx1, dy1),
                 c(dx2, dy2)))
  if (D==0) {
    return(c(Inf, Inf))
  }
  D1 <- det(rbind(P1, P2))
  D2 <- det(rbind(P3, P4))
  
  X <- det(rbind(c(D1, dx1),
                 c(D2, dx2)))/D
  Y <- det(rbind(c(D1, dy1),
                 c(D2, dy2)))/D
  
  if (interior.only) {
    ## Compute the fractions of L1 and L2 at which the intersection
    ## occurs
    lambda1 <- -((X-P1[1])*dx1 + (Y-P1[2])*dy1)/(dx1^2 + dy1^2)
    lambda2 <- -((X-P3[1])*dx2 + (Y-P3[2])*dy2)/(dx2^2 + dy2^2)
    if (!((lambda1>0) & (lambda1<1) &
          (lambda2>0) & (lambda2<1))) {
      return(c(NA, NA))
    }
  }
  return(c(X, Y))
}

##' This is simlar to unique(), but spares rows which are duplicated, but 
##' at different points in the matrix
##' 
##' @title Remove identical consecutive rows from a matrix
##' @param P Source matrix
##' @return Matrix with identical consecutive rows removed.
##' @author David Sterratt
remove.identical.consecutive.rows <- function(P) {
  for (i in 2:nrow(P)) {
    if (identical(P[i-1,], P[i,])) {
      return(remove.identical.consecutive.rows(P[-i,]))
    }
  }
  return(P)
}

##' Suppose segments AB and CD intersect.  Point B is replaced by the
##' intersection point, defined B'.  Point C is replaced by a point C'
##' on the line B'D. The maxium distance of B'C' is given by the
##' parameter d. If the distance l B'D is less than 2d, the distance
##' B'C' is l/2.
##'
##' @title Remove intersections between adjacent segements in a closed path
##' @param P The points, as a 2-column matrix
##' @param d Criterion for maximum distance when points are inser
##' @return A new closed path without intersections
##' @author David Sterratt
##' @export
remove.intersections <- function(P, d=50) {
  N <- nrow(P)
  for (i in 1:N) {
    R <- line.line.intersection(P[i,],            P[mod1(i+1, N),],
                                P[mod1(i+2, N),], P[mod1(i+3, N),],
                                interior.only=TRUE)
    if (identical(P[mod1(i+1, N),], P[mod1(i+2, N),])) {
      R <- P[mod1(i+1, N),]
    }
    if (is.finite(R[1])) {
      message("Intersection found. Old points:")
      message(paste(" ", (P[i,])))
      message(paste(" ", (P[mod1(i+1, N),])))
      message(paste(" ", (P[mod1(i+2, N),])))
      message(paste(" ", (P[mod1(i+3, N),])))

      P[mod1(i+1, N),] <- R
      message("Point i+1 has been changed:")
      message(paste(" ", (P[i,])))
      message(paste(" ", (P[mod1(i+1, N),])))
      message(paste(" ", (P[mod1(i+2, N),])))
      message(paste(" ", (P[mod1(i+3, N),])))

      l <- vecnorm(P[mod1(i+1, N),] - P[mod1(i+3, N),])
      if (l > 2*d) {
        a <- d/l
      } else {
        a <- 0.5
      }
      message(paste(" ", (paste("a=", a))))
      message(paste(" ", (paste("l=", l))))
      P[mod1(i+2, N),] <- a*P[mod1(i+1, N),] + (1-a)*P[mod1(i+3, N),]
      message("New points:")
      message(paste(" ", (P[i,])))
      message(paste(" ", (P[mod1(i+1, N),])))
      message(paste(" ", (P[mod1(i+2, N),])))
      message(paste(" ", (P[mod1(i+3, N),])))
    }
  }
  return(P)
}


##' Return points on the unit circle in an anti-clockwise
##' direction. If \code{L} is not specified \code{n} points are
##' returned. If \code{L} is specified, the same number of points are
##' returned as there are elements in \code{L}, the interval between
##' successive points being proportional to \code{L}.
##'
##' @title Return points on the unit circle
##' @param n Number of points
##' @param L Intervals between points
##' @return The cartesian coordinates of the points
##' @author David Sterratt
circle <- function(n=12, L=NULL) {
  if (is.null(L)) {
    angles <- (0:(n-1))/n*2*pi
  } else {
    angles <- (cumsum(L)-L[1])/sum(L)*2*pi
  }
  return(cbind(x=cos(angles), y=sin(angles)))
}

##' Find the interections of the plane defined by the normal \code{n} and the
##' distance \code{d} expressed as a fractional distance along the side of
##' each triangle.
##'
##' @title Find the intersection of a plane with edges of triangles on
##' a sphere
##' @param phi Lattitude of grid points on sphere centred on origin.
##' @param lambda Longitude of grid points on sphere centred on origin.
##' @param T Triangulation
##' @param n Normal of plane
##' @param d Distance of plane along normal from origin.
##' @return Matrix with same dimensions as \code{T}. Each row gives
##' the intersection of the plane  with the corresponding triangle in
##' \code{T}. Column 1 gives the fractional distance from vertex 2 to
##' vertex 3. Column 2 gives the fractional distance from vertex 3 to
##' vertex 1. Column 2 gives the fractional distance from vertex 1 to
##' vertex 2. A value of \code{NaN} indicates that the corresponding
##' edge lies in the plane. A value of \code{Inf} indicates that the
##' edge lies parallel to the plane but outside it.
##' @author David Sterratt
compute.intersections.sphere <- function(phi, lambda, T, n, d) {
  P <- cbind(cos(phi)*cos(lambda),
             cos(phi)*sin(lambda),
             sin(phi))
  return(cbind((d - P[T[,2],] %*% n)/((P[T[,3],] - P[T[,2],]) %*% n),
               (d - P[T[,3],] %*% n)/((P[T[,1],] - P[T[,3],]) %*% n),
               (d - P[T[,1],] %*% n)/((P[T[,2],] - P[T[,1],]) %*% n)))
}


##' Convert locations of points on sphere in spherical coordinates to
##' points in 3D cartesian space
##'
##' @title Convert from spherical to Cartesian coordinates
##' @param phi vector of lattitudes of N points
##' @param lambda vector of longitudes of N points
##' @param R radius of sphere 
##' @return An N-by-3 matrix in which each row is the cartesian (X, Y,
##' Z) coordinates of each point
##' @author David Sterratt
sphere.spherical.to.sphere.cart <- function(phi, lambda, R=1) {
  P <- cbind(R*cos(phi)*cos(lambda),
             R*cos(phi)*sin(lambda),
             R*sin(phi))
  colnames(P) <- c("X", "Y", "Z")
  return(P)
}

##' Given a triangular mesh on a sphere described by mesh locations
##' (\code{phi}, \code{lambda}), a radius \code{R} and a triangulation
##' \code{Tt}, determine the Cartesian coordinates of points \code{cb}
##' given in barycentric coordinates with respect to the mesh.
##'
##' @title Convert barycentric coordinates of points in mesh on sphere
##' to cartesian coordinates 
##' @param phi Lattitudes of mesh points
##' @param lambda Longitudes of mesh points
##' @param R Radius of sphere 
##' @param Tt Triagulation
##' @param cb Object returned by tsearch containing information on the
##' triangle in which a point occurs and the barycentric coordinates
##' within that triangle
##' @return An N-by-3 matrix of the Cartesian coordinates of the points
##' @author David Sterratt
##' @export
bary.to.sphere.cart <- function(phi, lambda, R, Tt, cb) {
  ## Initialise output
  cc <- matrix(NA, nrow(cb$p), 3)
  colnames(cc) <- c("X", "Y", "Z")

  ## If there are no points, exit
  if (nrow(cb$p) == 0) {
    return(cc)
  }

  ## Find 3D coordinates of mesh points
  P <- sphere.spherical.to.sphere.cart(phi, lambda, R)

  ## Now find locations cc of datapoints in Cartesian coordinates  
  for(i in 1:nrow(cb$p)) {
    cc[i,] <- bary2cart(P[Tt[cb$idx[i],],], cb$p[i,])
  }
  return(cc)
}

##' Convert locations on the surface of a sphere in cartesian
##' (X, Y, Z) coordinates to spherical (phi, lambda) coordinates. 
##'
##' It is assumed that all points are lying on the surface of a sphere
##' of radius R.
##' @title Convert from Cartesian to spherical coordinates
##' @param P locations of points on sphere as N-by-3 matrix with
##' labelled columns "X", "Y" and "Z"
##' @param R radius of sphere 
##' @return N-by-2 Matrix wtih columns ("phi" and "lambda") of
##' locations of points in spherical coordinates 
##' @author David Sterratt
##' @export
sphere.cart.to.sphere.spherical <- function(P, R=1) {
  return(cbind(phi   =asin(P[,"Z"]/R),
               lambda=atan2(P[,"Y"], P[,"X"])))
}

##' Project spherical coordinate system \eqn{(\phi, \lambda)} to a polar
##' coordinate system \eqn{(\rho, \lambda)} such that the area of each
##' small region is preserved.
##'
##' This requires \deqn{R^2\delta\phi\cos\phi\delta\lambda =
##' \rho\delta\rho\delta\lambda}.  Hence \deqn{R^2\int^{\phi}_{-\pi/2}
##' \cos\phi' d\phi' = \int_0^{\rho} \rho' d\rho'}.  Solving gives
##' \eqn{\rho^2/2=R^2(\sin\phi+1)} and hence
##' \deqn{\rho=R\sqrt{2(\sin\phi+1)}}.
##' 
##' As a check, consider that total area needs to be preserved.  If
##' \eqn{\rho_0} is maximum value of new variable then
##' \eqn{A=2\pi R^2(\sin(\phi_0)+1)=\pi\rho_0^2}. So
##' \eqn{\rho_0=R\sqrt{2(\sin\phi_0+1)}}, which agrees with the formula
##' above.
##' @title Convert lattitude on sphere to radial variable in
##' area-preserving projection
##' @param phi Lattitude
##' @param R Radius
##' @return Coordinate \code{rho} that has the dimensions of length
##' @author David Sterratt
spherical.to.polar.area <- function(phi, R=1) {
  return(R*sqrt(2*(1 + sin(phi))))
}

##' This is the inverse of \code{\link{polar.cart.to.sphere.spherical}}
##'
##' @title Convert spherical coordinates on sphere to  polar
##' projection in Cartesian coordinates
##' @param r 2-column Matrix of spherical coordinates of points on
##' sphere. Column names are \code{phi} and \code{lambda}.
##' @param pa If \code{TRUE}, make this an area-preserving projection
##' @param preserve Quantity to preserve locally in the
##' projection. Options are \code{lattitude}, \code{area} or
##' \code{angle}
##' @return 2-column Matrix of Cartesian coordinates of points on polar
##' projection. Column names should be \code{x} and \code{y}
##' @author David Sterratt
##' @export
sphere.spherical.to.polar.cart <- function(r, pa=FALSE, preserve="lattitude") {
  ## FIXME: This function should be deprecated in favour of
  ## azimuthal.equalarea and azimuthal.equidistant in projections.R
  rho <- NULL
  if (pa)
    preserve <- "area"
  if (preserve=="area") {
    rho <- sqrt(2*(1 + sin(r[,"phi"])))
    ## rho <- spherical.to.polar.area(r[,"phi"])
  }
  if (preserve=="angle") {
    ## rho = alpha*sqrt(2*(1+sin(phi))/(1-sin(phi)))
    rho <- sqrt(2*(1+sin(r[,"phi"]))/(1-sin(r[,"phi"])))
  }
  if (preserve=="lattitude") {
    rho <- pi/2 + r[,"phi"]
  }
  if (is.null(rho))
    stop(paste("preserve argument", preserve, "not recognised"))
  x <- rho*cos(r[,"lambda"])
  y <- rho*sin(r[,"lambda"])
  return(cbind(x=x, y=y))
}

##' This is the inverse of \code{\link{sphere.spherical.to.polar.cart}}
##'
##' @title Convert polar projection in Cartesian coordinates to
##' spherical coordinates on sphere
##' @param r 2-column Matrix of Cartesian coordinates of points on
##' polar projection. Column names should be \code{x} and \code{y}
##' @param pa If \code{TRUE}, make this an area-preserving projection
##' @param preserve Quantity to preserve locally in the
##' projection. Options are \code{lattitude}, \code{area} or
##' \code{angle}
##' @return 2-column Matrix of spherical coordinates of points on
##' sphere. Column names are \code{phi} and \code{lambda}.
##' @author David Sterratt
##' @export
polar.cart.to.sphere.spherical <- function(r, pa=FALSE, preserve="lattitude") {
  ## FIXME: This function should be deprecated in favour of as-yet
  ## unwritten functions azimuthal.equalarea.inverse and
  ## azimuthal.equidistant.inverse in projections.R
  rho <- NULL
  if (pa)
    preserve <- "area"
  rho2 <- r[,"x"]^2 + r[,"y"]^2
  if (preserve=="area") {
    ## Need to make sure that the argument is not greater that 1. This
    ## can happen when passing a values produced from a Cartesian grid
    sinphi <- rho2/2 - 1
    sinphi[sinphi>1] <- NA
    phi <- asin(sinphi)
  }
  if (preserve=="angle") {
    ## phi = asin((rho^2/alpha^2-2)/(rho^2/alphpa^2+2))
    phi <- asin((rho2 - 2)/(rho2 + 2))
  }
  if (preserve=="lattitude") {
    phi <- sqrt(rho2) - pi/2
  }
  if (is.null(phi))
    stop(paste("preserve argument", preserve, "not recognised"))
  lambda <- atan2(r[,"y"], r[,"x"])
  return(cbind(phi=phi, lambda=lambda))
}

##' @title Convert azimuth-elevation coordinates to spherical coordinates
##' @param r Coordinates of points in azimuth-elevation coordinates
##' represented as  2 column matrix with column names \code{alpha}
##' (elevation) and \code{theta} (azimuth).
##' @param r0 Direction of the axis of the sphere on which to project
##' represented as a 2 column matrix of with column names \code{alpha}
##' (elevation) and \code{theta} (azimuth).
##' @return 2-column matrix of spherical coordinates of points with
##' column names \code{psi} (colattidude) and \code{lambda} (longitude).
##' @author David Sterratt
##' @examples
##' r0 <- cbind(alpha=0, theta=0)
##' r <- rbind(r0, r0+c(1,0), r0-c(1,0), r0+c(0,1), r0-c(0,1))
##' azel.to.sphere.colattitude(r, r0)
##' @export
azel.to.sphere.colattitude <- function(r, r0) {
  ## Find Cartesian coordinates of aziumuth and elevation on sphere
  rc <- cbind(cos(r[,"alpha"])*sin(r[,"theta"]),
              cos(r[,"alpha"])*cos(r[,"theta"]),
              sin(r[,"alpha"]))

  ## Find Cartesian coordinates of aziumuth and elevation on sphere
  r0c <- cbind(cos(r0[,"alpha"])*sin(r0[,"theta"]),
               cos(r0[,"alpha"])*cos(r0[,"theta"]),
               sin(r0[,"alpha"]))

  ## Angle made with optic axis is
  psi <- acos(rc %*% t(r0c))

  ## Projection onto the plane perpendicular to the optic axis
  pc <- rbind(cbind(-cos(r0[,"theta"]),
                     sin(r0[,"theta"]),
                     0),
              cbind(-sin(r0[,"alpha"])*sin(r0[,"theta"]),
                    -sin(r0[,"alpha"])*cos(r0[,"theta"]),
                     cos(r0[,"alpha"]))) %*% t(rc)
  print(r0c)
  print(rc)
  print(t(pc))
  lambdap <- atan2(pc[2,], pc[1,])

  out <- cbind(psi, lambdap)
  colnames(out) <- c("psi", "lambda")
  
  return(out)
}

##' This rotates points on sphere by specifying the direction its
##' polar axis, i.e. the axis going through (90, 0), should point
##' after (a) a rotation about an axis through the poins (0, 0) and
##' (0, 180) and (b) rotation about the original polar axis.
##' @title Rotate axis of sphere
##' @param r Coordinates of points in spherical coordinates
##' represented as  2 column matrix with column names \code{phi}
##' (lattitude) and \code{lambda} (longitude).
##' @param r0 Direction of the polar axis of the sphere on which to project
##' represented as a 2 column matrix of with column names \code{phi}
##' (lattitude) and \code{lambda} (longitude).
##' @return 2-column matrix of spherical coordinates of points with
##' column names \code{phi} (lattidude) and \code{lambda} (longitude).
##' @author David Sterratt
##' @examples
##' r0 <- cbind(phi=0, lambda=-pi/2)
##' r <- rbind(r0, r0+c(1,0), r0-c(1,0), r0+c(0,1), r0-c(0,1))
##' r <- cbind(phi=pi/2, lambda=0)
##' rotate.axis(r, r0)
##' @export
rotate.axis <- function(r, r0) {
  ## Find cartesian coordinates of points on sphere
  P <- sphere.spherical.to.sphere.cart(r[,"phi"], r[,"lambda"])

  ## If we are not changing the longitude of the main axis, do nothing
  ## apart from convert back to spherical coordinates, which has the
  ## happy sideeffect of normalising the angles within the range (-pi,
  ## pi].
  if (r0[,"phi"] != pi/2) {
    ## Rotate them about the equatorial axis through the 0 degrees meridian
    ## (the x-axis)
    dp <- pi/2 - r0[,"phi"]
    P <- P %*% rbind(c(1, 0, 0),
                     c(0,  cos(dp), sin(dp)),
                     c(0, -sin(dp), cos(dp)))
    ## This will have taken the North pole to (0, -90). Hence we need to
    ## rotate by another 90 degrees to get to where we want to.
    
    ## Then rotate about the z-axis
    dl <- r0[,"lambda"] + pi/2
    P <- P %*% rbind(c( cos(dl), sin(dl), 0),
                     c(-sin(dl), cos(dl), 0),
                     c(0, 0, 1))
    colnames(P) <- c("X", "Y", "Z")
  }
  return(sphere.cart.to.sphere.spherical(P))
}

##' @title Bring angle into range
##' @param theta Angle to bring into range \code{[-pi, pi]}
##' @return Normalised angle
##' @author David Sterratt
##' @export
normalise.angle <- function(theta) {
  i <- which((theta < -pi) | (theta > pi) & !is.na(theta))
  if (length(i) > 0) {
    theta[i] <- ((theta[i] + pi) %% (2*pi)) - pi
  }
  return(theta)
}

##' This in the inverse of \code{\link{sphere.cart.to.sphere.wedge}}
##'
##' @title Convert from 'wedge' to Cartesian coordinates
##' @param psi vector of slice angles of N points
##' @param f vector of fractional distances of N points
##' @param phi0 rim angle as colatitude
##' @param R radius of sphere 
##' @return An N-by-3 matrix in which each row is the cartesian (X, Y,
##' Z) coordinates of each point
##' @export
##' @author David Sterratt
sphere.wedge.to.sphere.cart <- function(psi, f, phi0, R=1) {
  r <- sqrt(sin(phi0)^2 + cos(phi0)^2*cos(psi)^2)
  y0 <- -sin(psi)*cos(psi)*cos(phi0)
  z0 <- -sin(psi)*sin(psi)*cos(phi0)
  alpha0 <- asin(sin(phi0)/r)
  alpha <- alpha0 + f*(2*pi - 2*alpha0)
  P <- cbind(R*r*sin(alpha),
             R*(y0 - r*sin(psi)*cos(alpha)),
             R*(z0 + r*cos(psi)*cos(alpha)))
  colnames(P) <- c("X", "Y", "Z")
  return(P)
}

##' Convert points in 3D cartesian space to locations of points on
##' sphere in 'wedge' coordinates (\var{psi}, \var{f}).  Wedges are
##' defined by planes inclined at an angle \var{psi} running through a
##' line between poles on the rim above the x axis.  \var{f} is the
##' fractional distance along the circle defined by the intersection
##' of this plane and the curtailed sphere.
##'
##' @title Convert from Cartesian to 'wedge' coordinates
##' @param P locations of points on sphere as N-by-3 matrix with
##' labelled columns "X", "Y" and "Z"
##' @param phi0 rim angle as colatitude
##' @param R radius of sphere 
##' @return 2-column Matrix of 'wedge' coordinates of points on
##' sphere. Column names are \code{phi} and \code{lambda}.
##' @export
##' @author David Sterratt
sphere.cart.to.sphere.wedge <- function(P, phi0, R=1) {
  ## Test that points lie approximately on the unit sphere Note the
  ## coordinates produced by bary.to.sphere.cart do not lie exactly on
  ## a sphere; they lie on the triangles that approximate the sphere.
  Rs <- sqrt(rowSums(P^2))
  if (any(abs(Rs - R)/R > 0.1)) {
    print(abs(Rs - R)/R)
    stop("Points do not lie approximately on unit sphere")
  }
  ## Normalise to unit sphere
  P <- P/Rs
  ## Wedge angle, making sure this lies within [-pi/2, pi/2]
  psi <- atan2(P[,"Y"], - cos(phi0) - P[,"Z"])
  psi[psi >  pi/2 + 1E-6] <- psi[psi >  pi/2 + 1E-6] - pi
  psi[psi < -pi/2 - 1E-6] <- psi[psi < -pi/2 - 1E-6] + pi
  r <- sqrt(sin(phi0)^2 + cos(phi0)^2*cos(psi)^2)
  y0 <- -sin(psi)*cos(psi)*cos(phi0)
  z0 <- -sin(psi)*sin(psi)*cos(phi0)
  v <- -(P[,"Y"] - y0)*sin(psi) + (P[,"Z"] - z0)*cos(psi)
  alpha <- atan2(P[,"X"], v)
  ## Make sure angles in the second quadrant are negative
  ## FIXME: we could avoid this by defining the angle \alpha differently
  inds <- alpha<0
  alpha[inds] <- alpha[inds] + 2*pi
  alpha0 <- asin(sin(phi0)/r)
  f <- (alpha - alpha0)/(2*pi - 2*alpha0)
  # return(cbind(psi=psi, f=f, v=v, alpha0, alpha, y0))
  Pw <- cbind(psi=psi, f=f)
  rownames(Pw) <- NULL
  ## Check that f is in bounds
  if (any(f > 1) || any(f < 0)) {
    print(Pw)
    stop("f is out of bounds")
  }
  return(Pw)
}

##' Convert points in 3D cartesian space to locations of points on
##' sphere in 'dualwedge' coordinates (\var{fx}, \var{fy}).  Wedges
##' are defined by planes inclined at angle running through a line
##' between poles on the rim above the x axis or the y-axis.  \var{fx}
##' and \var{fy} are the fractional distances along the circle defined
##' by the intersection of this plane and the curtailed sphere.
##'
##' @title Convert from Cartesian to 'dualwedge' coordinates
##' @param P locations of points on sphere as N-by-3 matrix with
##' labelled columns "X", "Y" and "Z"
##' @param phi0 rim angle as colatitude
##' @param R radius of sphere 
##' @return 2-column Matrix of 'wedge' coordinates of points on
##' sphere. Column names are \code{phi} and \code{lambda}.
##' @export
##' @author David Sterratt
sphere.cart.to.sphere.dualwedge <- function(P, phi0, R=1) {
  Pwx <- sphere.cart.to.sphere.wedge(P, phi0, R=R)
  Pwy <- sphere.cart.to.sphere.wedge(cbind(X=P[,"Y"], Y=-P[,"X"], Z=P[,"Z"]),
                                     phi0, R=R)
  Pw <- cbind(fx=Pwx[,"f"], fy=Pwy[,"f"])
  rownames(Pw) <- NULL
  return(Pw)
}


##' On a sphere the central angle between two points is defined as the
##' angle whose vertex is the centre of the sphere and that subtends
##' the arc formed by the great circle between the points. This
##' function computes the central angle for two points \eqn{(\phi_1,
##' \lambda_1)}{(phi1, lambda1)} and \eqn{(\phi_2,\lambda_2)}{(phi2,
##' lambda2)}.
##'
##' @title Central angle between two points on a sphere
##' @param phi1 Lattitude of first point
##' @param lambda1 Longitude of first point
##' @param phi2 Lattitude of second point
##' @param lambda2 Longitude of secone point
##' @return Central angle
##' @source Wikipedia \url{http://en.wikipedia.org/wiki/Central_angle}
##' @author David Sterratt
##' @export
central.angle <- function(phi1, lambda1, phi2, lambda2) {
  return(acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2)))
}

##' The Karcher mean of a set of points on a manifold is defined as
##' the point whose sum of squared Riemmann distances to the points is
##' minimal. On a sphere using sphereical coordinates this distance
##' can be computed using the formula for central angle.
##'
##' @title Karcher mean on the sphere
##' @param x Matrix of points on sphere as N-by-2 matrix with labelled
##' columns \\code{phi} (lattitude) and \code{lambda} (longitude)
##' @param na.rm logical value indicating whether \code{NA} values should
##' be stripped before the computation proceeds.
##' @param var logical value indicating whether variance should be
##' returned too.
##' @return Vector of means with components named \code{phi} and
##' \code{lambda}. If \code{var} is \code{TRUE}, a list containing
##' mean and varience in elements \code{mean} and \code{var}.
##' @references Heo, G. and Small, C. G. (2006). Form representations
##' and means for landmarks: A survey and comparative
##' study. \emph{Computer Vision and Image Understanding},
##' 102:188-203.
##' @seealso \code{\link{central.angle}}
##' @author David Sterratt
##' @export
karcher.mean.sphere <- function(x, na.rm=FALSE, var=FALSE) {
  if (na.rm) {
    x <- na.omit(x)
  }
  ## Default output values, if there x has zero rows
  mu     <- c(phi=NA, lambda=NA)
  sigma2 <- c(phi=NA, lambda=NA)
  ## x has one row - needed to prevent crash
  if (nrow(x) == 1) {
    mu     <- c(x[1,])
    sigma2 <- c(phi=NA, lambda=NA)
  }
  ## x has more than one row
  if (nrow(x) >= 2) {
    ## Compute first estimate of mean by computing centroid in 3D and
    ## then finding angle to this
    P <- cbind(cos(x[,"phi"])*cos(x[,"lambda"]),
               cos(x[,"phi"])*sin(x[,"lambda"]),
               sin(x[,"phi"]))
    N <- nrow(P)
    P.mean <- apply(P, 2, mean)
    phi.mean <-    asin(P.mean[3])
    lambda.mean <- atan2(P.mean[2], P.mean[1])

    ## Now minimise sum of squared distances
    if (all(!is.nan(c(phi.mean, lambda.mean)))) {
      opt <- optim(c(phi.mean, lambda.mean),
                   function(p) { sum((central.angle(x[,"phi"], x[,"lambda"], p[1], p[2]))^2) })
      mu <- opt$par
      names(mu) <- c("phi", "lambda")
      sigma2 <- opt$value/N
    } else {
      mu <-     cbind(phi=NaN, lambda=NaN)
      sigma2 <- cbind(phi=NaN, lambda=NaN)
    }
  }
  ## Assemble output
  if (var) {
    X <- list(mean=mu, var=sigma2)
  } else {
    X <- mu
  }
  return(X)
}

##' @title Create grid on projection of hemisphere onto plane
##' @param pa If \code{TRUE}, make this an area-preserving projection
##' @param res Resolution of grid
##' @param phi0 Value of \code{phi0} at edge of grid
##' @return List containing:
##' \item{\code{s}}{Grid locations in spherical coordinates}
##' \item{\code{c}}{Grid locations in Cartesian coordinates on plane}
##' \item{\code{xs}}{X grid line locations in Cartesian coordinates on plane}
##' \item{\code{ys}}{Y grid line locations in Cartesian coordinates on plane}
##' @author David Sterratt
##' @export
create.polar.cart.grid <- function(pa, res, phi0) {
  lim <- sphere.spherical.to.polar.cart(cbind(phi=phi0, lambda=0), pa)[1,"x"]
  xs <- seq(-lim, lim, len=res)
  ys <- seq(-lim, lim, len=res)

  ## Create grid
  gxs <- outer(xs, ys*0, "+")
  gys <- outer(xs*0, ys, "+")

  ## gxs and gys are both res-by-res matrices We now combine both
  ## matrices as a res*res by 2 matrix. The conversion as.vector()
  ## goes down the columns of the matrices gxs and gys
  gc <- cbind(x=as.vector(gxs), y=as.vector(gys))

  ## Now convert the cartesian coordinates to polar coordinates
  gs <- polar.cart.to.sphere.spherical(gc, pa)
  return(list(s=gs, c=gc, xs=xs, ys=ys))
}

