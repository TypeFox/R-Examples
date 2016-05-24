# Package SphericalCubature
# John P. Nolan, April 2013, jpnolan@american.edu
# This package computes integrals over spheres, balls, and
#    annular regions in n-dimensions.
#
# Functions integratePolynomialSphere and integratePolynomialBall
#     exactly integrate polynomials over a sphere/ball.
# Function integrateSphereStroud11 uses Stroud 11 degree rule
#     to integrate over the unit sphere. It is faster but only
#     reliable for smooth functions.
# Functions adaptIntegrateSphere and adaptIntegrateBall use
#     adaptive integration using the multivariate polar representation
#     and the R package cubature, which performs multivariate adaptive
#     integration over that rectangular region.  Both functions have
#     extensions, named adaptIntegrateSphereSplit and adaptIntegrateBallSplit,
#     which allow you to manually split the region of integration to
#     capture spikes, etc. in the integrand.
# Utility functions:
#     rect2polar and polar2rect convert to/from polar coordinates in
#         n-dimensions
#     nextGraySubset and nextMultiIndex calculate particular subsets
#         needed in the computations.
#
########################################################################
#  n dimensional polar coordinates are given by the following:
#  rectangular x=(x[1],...,x[n])  corresponds to
#     polar  (r,phi[1],...,phi[n-1]) by
#  x[1]  = r*cos(phi[1])
#  x[2]  = r*sin(phi[1])*cos(phi[2])
#  x[3]  = r*sin(phi[1])*sin(phi[2])*cos(phi[3])
#   ...
#  x[n-1]= r*sin(phi[1])*sin(phi[2])*...*sin(phi[n-2])*cos(phi[n-1])
#  x[n]  = r*sin(phi[1])*sin(phi[2])*...*sin(phi[n-2])*sin(phi[n-1])
#  Here phi[1],...,phi[n-2] in [0,pi), and phi[n-1] in [0,2*pi).
#  For multivariate integration, the Jacobian of the above tranformation
#    is J(phi) = r^(n-1) * prod( sin(phi[1:(n-2)])^((n-2):1) );
#    note that phi[n-1] does not appear in the Jacobian.
########################################################################
require(cubature) # for adaptIntegrate*( ) functions
########################################################################
integrateSpherePolynomial <- function( p, valueOnly=TRUE ) {
# compute the exact integral over the sphere in n dimensions of
# a polynomial  p(x[1],...,x[n])=sum (coef[i] * x[1]^k[i,1] * ... * x[n]^k[i,n]),
# where the sum is over i=1,...,m.  This is specified as a list
# with fields
#    p$coef[1:m] vector of doubles
#    p$k[1:m,1:n] matrix of integers
#    m and n a given implicitly in the sizes of these arrays
# output - a list with two fields:
#    $value - a double containing the value of the integral
#    $term - a vector of length m of values used in the computation
#          and which are used by function IntegratePolynomialBall( ... )
#
# Method is from Folland (2001), MAA Monthly 108, pg. 446-448.

m <- length(p$coef)
storage.mode(p$k) <- "integer"
if (m != nrow(p$k)) {
  stop("Error in integratePolynomialSphere: length(p$coef) should be the same as he number of rows of p$k\n")
}
if ( any(p$k < 0) ) { stop("Error in integratePolynomialSphere: some p$k[i,j] is negative") }
n <- ncol(p$k)

# loop through the terms in the polynomial, integrating each monomial exactly
term <- double(m)
for (i in 1:m) {
  ki <- p$k[i, ]
  if ( !any( ki %% 2L )) { # test that all powers are even, otherwise term=0
    beta <- (ki+1.0)/2.0
    term[i] <- 2*p$coef[i]*prod(gamma(beta))/gamma(sum(beta))
  }
}
value <- sum(term)
if (valueOnly) return(value)
else return( list(value=value,term=term) )}
########################################################################
integrateBallPolynomial <- function( p, R=c(0.0,1.0) ) {
# compute the exact integral over (part of) a ball in n dimensions of
# a polynomial p(x[1],...,x[n]), which is specified as in the
# function IntegratePolynomialSphere( ).  The default is to integrate
# over the unit ball:  0 <= r <= 1, but one can specify inner and outer
# r bounds to get integration over the annulus: R[1] <= r <= R[2].
# Method is from reference in IntegratePolynomialSphere( ).

n <- ncol(p$k)
b <- integrateSpherePolynomial( p, valueOnly=FALSE )
value <- 0.0
for (i in 1:length(p$coef)) {
  power <- n + sum(p$k[i,])
  value <- value + b$term[i]*(R[2]^power-R[1]^power)/power
}
return(value)}
########################################################################
adaptIntegrateSphere <- function( f, n, lowerLimit=rep(0.0,n-1),
                    upperLimit=c(rep(pi,n-2),2*pi), tol=1e-05, ... ) {
# Integrate f(x[1],...,x[n]) over the portion of the unit sphere
#     bounded by *polar* coordinates lowerLimit and upperLimit.
#     If the integrand has abrupt changes, the adaptive integration
#     routine may miss that change.  You can prevent this if you
#     know the location of the change points by calling
#     adaptIntegrateSphereSplit( ), which will split up the region
#     of integration at/near specified points.
#
# The result is a list with multiple fields.  In all cases, the value of
# the integral is stored in the field $value.  Other fields are specific to
# the dimension: if n=2, the list is what is returned by integrate( );
# if n > 2, the list is what is returned by adaptIntegrate( ) from
# the package cubature.

n <- as.integer(n)
a <- adaptIntegrateCheck( n, lowerLimit, upperLimit, R=c(0.0,1.0), xstar=matrix(0.0,nrow=n,ncol=0), width=0.0 )

if (n == 2) {   # use adaptive univariate integration
  polar.f <- function( phi, rect.f, ... ) {
      # function integrate( ) in base R is vectorized; rect.f( ) may not
      #   be vectorized, so we evaluate f on argument at a time.  This will
      #   be slower, but we do it to be consistent with n > 2.
      y <- double(length(phi))
      for (i in 1:length(phi)) {
        x <- polar2rect( 1.0, phi[i])
        y[i] <- rect.f(x,...)
      }
      return(y) }
    
  z <- integrate( polar.f, lower=lowerLimit, upper=upperLimit, rel.tol=tol, rect.f=f,... )
  if (z$message != "OK") warning(paste("Error in adaptIntegrateSphere: error return from integrate(...) ",z$message))
  return( z )
}

# From now on, dimension n > 2, use adaptIntegrate from package cubature
# internal integrand function
polar.f <- function( phi, rect.f, ... ) {
  x <- polar2rect( 1.0, phi )
  y <- rect.f(x,...)  # evaluate integrand
  
  # compute Jacobian from change to polar coordinates
  n <- length(phi)+1
  sint <- sin(phi)
  if (n > 2) { y <- y*prod( sint[1:(n-2)]^((n-2):1) ) }
  return( y ) }

z <- adaptIntegrate( polar.f,  lowerLimit=lowerLimit, upperLimit=upperLimit,
         tol=tol, rect.f=f, ... )
if(z$returnCode != 0) warning(paste("Error in adaptIntegrateSphere: return from adaptIntegrate=",z$returnCode))
z$value <- z$integral
return(z) }
########################################################################
adaptIntegrateSphereSplit <- function( f, n, xstar, width=0,
       lowerLimit=rep(0.0,n-1), upperLimit=c(rep(pi,n-2),2*pi), tol=1e-05, ... ) {
# Integrate f(x[1],...,x[n]) over the portion of the unit sphere
#     bounded by *polar* coordinates lowerLimit and upperLimit.
# If it is known that f(.) has abrupt changes (spike, cusp, etc).
#     in directions given by xstar[ ,j], j=1:m, then we can assist the
#     adaptive integration routine by splitting up the region of
#     integration into pieces, with rectangles of side width[j] around
#     each point xstar[,j].  (If width[j]=0, then just split at each
#     point xstar[ ,j].
# Since this function makes repetitive calls to adaptIntegrateSphere,
#     and that function prints a warning when something isn't right,
#     this function only returns a single number, the estimated
#     value of the integral

# check input and set up things
if (!is.matrix(xstar)) xstar <- as.matrix(xstar,ncol=1)
a <- adaptIntegrateCheck( n,  lowerLimit, upperLimit, R=c(0.0,1.0), xstar, width  )
m <- a$m
newWidth <- a$width

# partition the region of integration into subrectangles based on directions
# of the vectors xstar[ ,j], j=1:m
partition <- partitionRegion( xstar, newWidth, lowerLimit, upperLimit )

#for (i in 1:length(partition$count)) { cat(i,partition$count[i],"   ",partition$phi[[i]],"\n") }

value <- 0.0
size <- partition$count - 1
j <- rep(1L,n-1)
low <- rep(0.0,n-1)
up <- rep(0.0,n-1)
repeat {
  for (i in 1:(n-1)) {
    low[i] <- partition$phi[[i]][j[i]]
    up[i] <- partition$phi[[i]][j[i]+1]
  }
  #cat("j=",j,"    ",low,"   ",up,"\n")
  b <- adaptIntegrateSphere( f, n, low,up, tol, ... )
  value <- value + b$value
  j <- nextMultiIndex(j, size)
  if (j[1] < 0 ) { break }
}
return(value) }
########################################################################
adaptIntegrateBall <- function( f, n, lowerLimit=rep(0.0,n-1),
                    upperLimit=c(rep(pi,n-2),2*pi), R=c(0.0,1.0), tol=1e-05, ... ) {
# integrate f(x[1],...,x[n]) over the portion of the unit ball
#     bounded by *polar* coordinates lowerLimit and upperLimit and
#     radius R[1] <= r <= R[2].
#     If the integrand has abrupt changes, the adaptive integration
#     routine may miss that change.  You can prevent this if you
#     know the location of the change points by calling
#     adaptIntegrateBallSplit( ), which will split up the region
#     of integration at/near specified points.
#
# The result is a list with multiple fields.  In all cases, the value of
# the integral is stored in the field $value.  The list is what is
# returned by adaptIntegrate( ) from the package cubature.

n <- as.integer(n)
a <- adaptIntegrateCheck( n, lowerLimit, upperLimit, R=R, xstar=matrix(0.0,nrow=n,ncol=0), width=0.0 )

# use adaptIntegrate from package cubature
# internal integrand function
polar.f <- function( r.phi, rect.f, ... ) {
  n <- length(r.phi)
  x <- polar2rect( r.phi[1],r.phi[2:n] )
  y <- rect.f(x,...)  # evaluate integrand
  
  # compute Jacobian from change to polar coordinates

  sint <- sin(r.phi[2:n])
  if (n == 2L) { 
    y <- r.phi[1]*y
  } else { 
    y <- r.phi[1]^(n-1)*y*prod( sint[1:(n-2)]^((n-2):1) ) 
  }
  return( y ) }

z <- adaptIntegrate( polar.f,  lowerLimit=c(R[1],lowerLimit), upperLimit=c(R[2],upperLimit),
         tol=tol, rect.f=f, ... )
if(z$returnCode != 0) warning(paste("Error in adaptIntegrateBall: return from adaptIntegrate=",z$returnCode))
z$value <- z$integral
return(z) }
########################################################################
adaptIntegrateBallSplit <- function( f, n, xstar, width=0.0, lowerLimit=rep(0.0,n-1),
                    upperLimit=c(rep(pi,n-2),2*pi), R=c(0.0,1.0), tol=1e-05, ... ) {
# integrate f(x[1],...,x[n]) over the portion of the unit ball
#    bounded by *polar* coordinates lowerLimit and upperLimit and
#    radius R[1] <= r <= R[2]
# If it is known that f(.) has abrupt changes (spike, cusp, etc).
#     in directions given by xstar[ ,j], j=1:m, then we can assist the
#     adaptive integration routine by splitting up the region of
#     integration into pieces, with rectangles of side width[j] around
#     each point xstar[,j].  (If width[j]=0, then just split at each
#     point xstar[ ,j].
# Since this function makes repetitive calls to adaptIntegrateSphere,
#     and that function prints a warning when something isn't right,
#     this function only returns a single number, the estimated
#     value of the integral

# check input and set up things
if (!is.matrix(xstar)) xstar <- as.matrix(xstar,ncol=1)
a <- adaptIntegrateCheck( n,  lowerLimit, upperLimit, R=R, xstar=xstar, width=width  )
m <- a$m
newWidth <- a$width

# partition the region of integration into subrectangles based on directions
# of the vectors xstar[ ,j], j=1:m
partition <- partitionRegion( xstar, newWidth, lowerLimit, upperLimit )

#for (i in 1:length(partition$count)) { cat(i,partition$count[i],"   ",partition$phi[[i]],"\n") }

value <- 0.0
size <- partition$count - 1
j <- rep(1L,n-1)
low <- rep(0.0,n-1)
up <- rep(0.0,n-1)
repeat {
  for (i in 1:(n-1)) {
    low[i] <- partition$phi[[i]][j[i]]
    up[i] <- partition$phi[[i]][j[i]+1]
  }
  #cat("j=",j,"    ",low,"   ",up,"\n")
  b <- adaptIntegrateBall( f, n, low,up, R=R, tol, ... )
  value <- value + b$value
  j <- nextMultiIndex(j, size)
  if (j[1] < 0 ) { break }
}
return(value) }
########################################################################
polar2rect <- function( r, phi ) {
# Convert polar coordinates to rectangular coordinates in n-dimensions.
# If r is a scalar, then the polar point (r,phi[1:(n-1)]) is converted to
#     rectangular coordinates x[1:n].
# If r is a vector of length m, then phi should be a matrix of dimensions (n-1) x m,
#     and the result is a matrix x[1:n,1:m], with columns of x being the x
#     coordinates of points (r[j],phi[,j]).
# The result is always a matrix x of size (n x m).
#
m <- length(r)
if (!is.matrix(phi)) { phi <- as.matrix(phi,ncol=1) }
stopifnot( m == ncol(phi))
n <- nrow(phi) + 1
x <- matrix(0.0,nrow=n,ncol=m)
for (j in 1:m) {
  col.cos <- cos(phi[,j])
  col.sin <- sin(phi[,j])
  s <- c( col.cos[1], rep(col.sin[1],n-1) )
  if (n > 2) {
    for (k in 2:(n-1)) {
      s[k] <- s[k]*col.cos[k]
      s[(k+1):n] <- s[(k+1):n]*col.sin[k]
    }
  }
  x[,j] <- r[j]*s
}
return(x) }
########################################################################
rect2polar <- function( x ) {
# Convert from rectangular to polar coordinates in n dimensions.
# If x[1:n] is a vector in rectangular coordinates, convert to
#     polar coordinates (r,phi[1:(n-1)])
# If x[1:n,1:m] is a matrix, convert each column x[ ,i] to polar coordinates
#   given by r[i] and phi[,i]

if(!is.matrix(x)) { x <- as.matrix(x,ncol=1) }
n <- nrow(x)
m <- ncol(x)
r <- rep(0.0,m)
phi <- matrix(0.0,nrow=n-1,ncol=m)
for (j in 1:m) {
  rsq <- x[,j]^2
  cum.rsq <- cumsum(rev(rsq))
  r[j] <- sqrt( cum.rsq[n] )
  if (r[j] > 0.0) {
    if (n>2) {
      for (k in 1:(n-2)) {
        phi[k,j] <- atan2( sqrt(cum.rsq[n-k]), x[k,j] )
      }
    }
    phi[n-1,j] <- 2*atan2( x[n,j], x[n-1,j]+sqrt(cum.rsq[2] ) )
  }
}
return(list(r=r,phi=phi))}
########################################################################
########################################################################
#  Stroud integration and related functions, adapted from fortran and
#  C code by John Burkhart found at
#    http://people.sc.fsu.edu/~jburkardt/f77_src/stroud/stroud.html
#    http://people.sc.fsu.edu/~jburkardt/c_src/stroud/stroud.html
#  Based on the book by A. H. Stroud, Approximate Calculation of
#  multiple integrals, 1971, page 296-297.
########################################################################
########################################################################
integrateSphereStroud11 <- function( f, n, ... ) {
# numerically integrate the function f( x, ... ) over the unit sphere in R^n
# using an 11th degree Stroud cubature
#
# ... are arguments passed to f( )
#
# Based on routine sphere_unit_11_nd by Burkhardt.
# Changes:
#  (1) the outer loops for S31 and S32 in Fortran should be DO I=1,N-1.
#  (2) in dimension 7, Burkhardt's program a mistake
#     in constant coef21[7], it should be 0.0337329118818, not 0.0336329118818

# check that dimesion isn't too small or too large
if (n < 3) stop("Error in Stroud11IntegrateSphere: n < 3")
if (n > 16) stop("Error in Stroud11IntegrateSphere: n > 16")

# define quadrature coefficients
coef1 <- c( 0.0, 0.0, 0.128571428571, 0.0518518518518, 0.0211979378646, 0.281250000000,
            1.11934731935, 2.82751322751, 5.68266145619, 9.93785824515,  15.8196616478,
            23.5285714285, 33.2409299392, 45.1113811729, 59.2754264177, 75.8518518518)
coef21 <- c( 0.0, 0.0, 0.163795782462, 0.0967270533860, 0.0638253880175, 0.0452340041459,
             0.0337329118818, 0.0261275095270, 0.0208331595340, 0.0169937111647,
             0.0141147212492, 0.0118949128383, 0.0101424250926, 0.00873046796644,
             0.00757257014768, 0.00660813369775)
coef22 <- c( 0.0, 0.0, 0.126680408014, 0.0514210947621, 0.0213579471658, -0.108726067638,
            -0.371589499738, -0.786048144448, -1.36034060198, -2.09547695631, -2.98784764467,
            -4.03107480702, -5.21726499521, -6.53783099707, -7.98401677102, -9.54722261180)
coef31 <- c( 0.0, 0.0, 0.0, 0.0592592592592, 0.0236639091329, 0.0525940190875, 0.0925052768546,
             0.141316953438, 0.196818580052, 0.257027634179, 0.320299222258, 0.385326226441,
             0.451098131789, 0.516849445559, 0.582010515746, 0.646165210110)
coef32 <- c( 0.0, 0.0, 0.0, 0.0, 0.0316246294890, 0.0207194729760, 0.0144303800811,
             0.0105348984135, 0.00798435122193, 0.00623845929545, 0.00499896882962,
             0.00409176297655, 0.00341037426698, 0.00288710646943, 0.00247745182907,
             0.00215128820597)

quad <- 0.0
x <- rep(1.0/sqrt(n),n)  # initial point on unit sphere
  # S1
  gray <- list(more=FALSE,n=n)
  repeat {
    gray <- nextGraySubset( gray )
    if (gray$iadd != 0L) { x[gray$iadd] <- -x[gray$iadd] }
    quad <- quad + coef1[n] * f ( x, ... )
    if (!gray$more) {break}
  }

  # S21
  r1 <- ( ( n + 6 ) - 4.0 * sqrt ( 3.0 ) )  / ( n * n + 12 * n - 12 )
  r1 <- sqrt ( r1 )

  s1 <- ( ( 7 * n - 6 ) + ( 4 * ( n - 1 ) ) * sqrt ( 3.0 ) ) / ( n * n + 12 * n - 12 )
  s1 <- sqrt ( s1 )

  for( i in 1:n ) {
    x <- rep(r1,n)
    x[i] <- s1
    gray$more <- FALSE
    repeat {
      gray <-  nextGraySubset (gray )
      if ( gray$iadd != 0L ) { x[gray$iadd] = -x[gray$iadd] }
      quad <- quad + coef21[n] * f ( x, ... )
      if (!gray$more) {break}
    }
  }

  #  S22
  r2 <- ( ( n + 6 ) + 4.0 * sqrt ( 3.0 ) )  /  ( n * n + 12 * n - 12 )
  r2 <- sqrt ( r2 )
  s2 = ( ( 7 * n - 6 ) - ( 4 * ( n - 1 ) ) * sqrt ( 3.0 ) ) / ( n * n + 12 * n - 12 )
  s2 <- sqrt ( s2 )
  for (i in 1:n) {
    x <- rep(r2, n)
    x[i] <- s2
    gray$more <- FALSE
    repeat {
      gray <- nextGraySubset ( gray )
      if ( gray$iadd != 0L ) { x[gray$iadd] = -x[gray$iadd] }
      quad = quad + coef22[n] * f ( x, ... )
      if ( !gray$more ) {break}
    }
  }
  
  #  S31
  u1 <- ( ( n + 12 ) + 8.0 * sqrt ( 3.0 ) ) /  ( n * n + 24 * n - 48 )
  u1 <- sqrt ( u1 )
  v1 <- ( ( 7 * n - 12 ) - ( 4 * n - 8 ) * sqrt ( 3.0 ) ) / ( n * n + 24 * n - 48 )
  v1 <- sqrt ( v1 )
  for (i in 1:(n-1)) {  # change from fortran code
    for (j in (i+1):n) {
      x <- rep(u1,n)
      x[i] = v1
      x[j] = v1
      gray$more <- FALSE
      repeat {
        gray <- nextGraySubset ( gray )
        if ( gray$iadd != 0L ) { x[gray$iadd] = -x[gray$iadd] }
        quad <- quad + coef31[n]* f ( x, ... )
        if ( !gray$more ) { break }
      }
    }
  }

  #  S32
  u2 <- ( ( n + 12 ) - 8.0 * sqrt ( 3.0 ) ) /  ( n * n + 24 * n - 48 )
  u2 <- sqrt ( u2 )
  v2 <- ( ( 7 * n - 12 ) + ( 4 * n - 8 ) * sqrt ( 3.0 ) ) / ( n * n + 24 * n - 48 )
  v2 <- sqrt ( v2 )
  for (i in 1:(n-1)) {  # change from fortran code
    for (j in (i+1):n) {
      x <- rep(u2,n)
      x[i] <- v2
      x[j] <- v2
      gray$more <- FALSE
      repeat {
        gray <- nextGraySubset ( gray )
        if ( gray$iadd != 0L ) { x[gray$iadd] = -x[gray$iadd] }
        quad <- quad + coef32[n] * f ( x, ... )
        if ( !gray$more ) { break }
      }
    }
  }
  result <- quad * sphereArea(n) / 2.0^n
return( result) }
#############################################################################
sphereArea <- function( n, R=1.0 ) {
# calculate the surface area of the sphere {|x|=R} in R^n

area <- R^(n-1)*2.0*pi^(n/2)/gamma(n/2)

return(area) }
#############################################################################
ballVolume <- function( n, R=1.0 ) {
# calculate the volume of the ball {|x|<R} in n-dimensions

volume <- R^n * pi^(n/2)/gamma((n/2.0)+1.0)
return(volume)}
#############################################################################
#############################################################################
# remaining functions are internal, not meant to be called by user.
#############################################################################
#############################################################################
nextGraySubset <- function( gray.list ) {
# compute the next Gray subset
# The first time through, gray.list should be a list with two elements:
#    gray.list$more - initially set to FALSE
#    gray.list$n - always set to size of the set
# Returns a list gray.list which should be used in the next call to
#    this function (or you will not get the next Gray subset).
#    Main result is in gray.list$a, which is a sequence of 0s and 1s, indicating
#    which elements are in the set and which are not.  Other fields are
#    needed for successive calls to this function.

#   First set returned is the empty set.
new <- gray.list
if ( !new$more ) {
  new$iadd <- 1L
  new$ncard <- 0L
  new$more <- TRUE
  a <- rep(0L,new$n)
}  else {
  new$iadd <- 1L
  a <- new$a

  if ( ( new$ncard %% 2L ) != 0L ) {
    repeat {
      new$iadd <- new$iadd + 1L
      if ( a[new$iadd-1] != 0L ) { break;}
    }
  }

  a[new$iadd] <- 1 - a[new$iadd];
  new$ncard <- new$ncard + 2 * a[new$iadd] - 1;

  #  Last set returned is the singleton a[n]
  if ( new$ncard == a[new$n] ) { new$more <- FALSE }
}

new$a <- a
return(new) }
############################################################################
nextMultiIndex <- function( j, size ) {
# Find the next (in lexographical order) multi-index (j[1],...,j[n])
#  by stepping through all values:
#        j[1] in 1:size[1], j[2] in 1:size[2], ..., j[n] in 1:size[n]
# Caller should initialize j to all 1s: j <- rep(1L,n) before calling
#   this function for the first time, then successively pass the multi-index
#   returned by this function back to this function until all multi-indices
#   have been generated.  When all possible values have been enumerated,
#   j[1] is set to -1 and return.

n <- length(j)
stopifnot( n == length(size) )

k <- n  # start at last digit and increment
repeat {
  if (j[k] < size[k]) {
    # increment position k and return
    j[k] <- j[k]+1L
    break
  } else {
    if (k <= 1) {  # no more positions
      j[1] <- -1
      break
    }
    # no more choices in position k, reset that position to 1
    # and go to earlier position
    j[k] <- 1L
    k <- k - 1L
  }
}
return(j) }
########################################################################
partitionRegion <- function( xstar, width, lowerLimit, upperLimit ){
# partition a polar region by splitting the region around points given
#   by xstar[ ,j], j=1:m
# Output is a list with fields
#   $count[1:k] number of bounds in $phi[[i]],
#   $phi[[i]] angle bounds for i-th polar coordinate

# convert points in xstar to polar coordinates
a <- rect2polar( xstar )  # ignore a$r = lengths of vectors

n <- nrow(a$phi) + 1
m <- ncol(a$phi)
partition <- list(count=integer(n-1),phi=vector("list",0) )
for (i in 1:(n-1)) {
  newphi <- a$phi[i, ]
  newphi <- c(newphi-width,newphi+width)
  maxphi <- ifelse(i<n, pi, 2*pi)
  if (any(newphi < 0.0) ){
    j <- which(newphi < 0.0)
    newphi[j] <- newphi[j]+maxphi
  }
  if (any(newphi > maxphi)) {
    j <- which( newphi > maxphi )
    newphi[j] <- newphi[j]-maxphi
  }
  newphi <- c(newphi,lowerLimit[i],upperLimit[i])
  newphi <- pmax(newphi,lowerLimit[i])
  newphi <- pmin(newphi,upperLimit[i])
  newphi <- unique(sort(newphi))
  partition$count[i] <- length(newphi)
  partition$phi[[i]] <- newphi
}
return( partition ) }
########################################################################
adaptIntegrateCheck <- function(  n, lowerLimit, upperLimit, R, xstar, width  ) {
# check input arguments to adaptIntegrate*( ) functions.

n <- as.integer(n)
if (n < 2) stop("Error in adaptIntegrateCheck: dimension n < 2")
m <- ncol(xstar)
if (n != nrow(xstar) ) stop("Error in adaptIntegrateCheck: n is not equal nrow(xstar)")
if( length(lowerLimit) != (n-1))
  stop(paste("Error in adaptIntegrateCheck: length(lowerLimit) should be ",n-1) )
if( length(lowerLimit) != length(upperLimit) )
  stop("Error in adaptIntegrateCheck: length(lowerLimit != length(upperLimit)")
if (any(lowerLimit < 0)) stop("Error in adaptIntegrateCheck: lowerLimit has negative values")
if (any(lowerLimit > upperLimit)) stop("Error in adaptIntegrateCheck: some lowerLimit > upperLimit")
maxLimit <- c(rep(pi,n-2),2*pi)
if (any(upperLimit > maxLimit)) stop("Error in adaptIntegrateCheck: some upperLimit too large")

if (length(width) == 1L) { width <- rep(width,m) }
if (length(width) != m ) stop("Error in adaptIntegrateCheck: length(width) != m")
if (any(width > pi/2)) {
   warning("adaptIntegrateCheck: some width > pi/2; truncated to pi/2")
   j <- which(width > pi/2)
   width[j] <- pi/2
}
if (length(R) != 2) stop("adaptIntegrateCheck: length(R) must be 2")
if (R[1] < 0.0) stop("adaptIntegrateCheck: R[1] must be nonnegative")
if (R[1] >= R[2]) stop("adaptIntegrateCheck: R[2] must be bigger than R[1]")

return(list(m=m,width=width) )}
#######################################################################
