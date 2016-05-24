setupInterp <- function(x, y, doPoly=TRUE ) {
## Interpolation points are given as vectors x and y.
##    Preparation of Newton's (doPoly=TRUE)
##    or rational (doPoly=FALSE) Interpolation.
## Rational Interpolation is more flexible in
##  approximating functions with poles and asymptotes
##  Change of vectors [0..lengthX] to [1..lengthX+1]
  lengthX <- length( x )
  if (lengthX != length(y)) {
    stop(paste("vectors have different lengths: x:",lengthX,", y:",length(y)))
  }
  if ( doPoly )  { # mode <- Newton
    ii <- 1;  jj <- 1;  kk <- lengthX
    X <- Q <- T <- rep(NA,lengthX)
    while ( jj < kk ) { # interleave from below and above
      X[ii] <- x[jj];  Q[ii] <- y[jj];  ii <- ii+1;  jj <- jj+1
      X[ii] <- x[kk];  Q[ii] <- y[kk];  ii <- ii+1;  kk <- kk-1
    }
        ## uneven lengthX
    if ( jj == kk ) { X[ii] <- x[kk];  Q[ii] <- y[kk] }
    for ( ii in 1:lengthX)  { 
      y <- X[ii];  T[ii] <- Q[ii]
      for ( kk in seqm(ii-1,1,-1) ) {
        T[kk] <- (T[kk+1] - T[kk]) / (y - X[kk] )
      }
      Q[ii] <- T[1]
    }
  } else { # mode <- Thiele
    X <- x;  Q <- y
    for ( ii in 2:lengthX )  {## inverse differences *)
      yy <- X[ii-1];	q <- Q[ii-1]
     for (kk in (ii:lengthX))  {
  # anything/0 -> Inf is harmless at later stages, since 1/Inf = 0
    	Q[kk] <- ( yy - X[kk]) / (q - Q[kk])
     }
    }
  }
  return(list(x=X,q=Q,ModeI=doPoly))
} ## end  setupInterpolation

evalInterp <- function( xi, ss ) {
##* Give interpolated value at xi, using setup in ss from SetupInterpol
  lengthX <- length(ss$x)
  ff <- ss$q[lengthX]
  if ( ss$ModeI )  { ##  Newton
    for (kk in seqm(lengthX-1, 1, -1)) {
      ff <- ( xi-ss$x[kk] )*ff + ss$q[kk]
    }
  } else { ##  Thiele
    for (kk in seqm(lengthX-1, 1, -1)) {
      ff <- ( xi-ss$x[kk] )/ff + ss$q[kk]
    }
  }
  return ( ff )   
}  ## end  evalIntp

minInterp <- function( x, y, add = FALSE, doPoly=TRUE ) {
  # Interpolate (doPoly=TRUE .= Newton) three points
  # or four points (x,y) (Thiele),
  # For add=TRUE one more point is used.
  # Returns the argument of the numerical minimum if it exists,
  #   or NA if too few points are given.
  ks <- sort.list(x)
  x  <- x[ks];  f <- y[ks] 
  nn <- length(x)
  km <- which.min(f)
  bm <- 4+add-doPoly
  if (nn < bm) {
    if (add && nn == bm-1) { # add == FALSE would also work
      add <- FALSE
      bm  <- bm-1
    } else {
      return( list(xmin=NA, int=list(fct=NA, coef=NA, xx=NA)))
    }
  }
  if (nn >= bm) {
    k <- max(1,min(km-round(bm/2),nn-bm+1))
    x <- x[k:(k+bm-1)];  f <- f[k:(k+bm-1)]        
  }
  s <- setupInterp( x, f, doPoly )
  if ( doPoly )  { #  Newton
    if (add==0)  s$q[4] <- 0.0  # not needed: s$x[3] <- 0.0
    a <- s$x[1] + s$x[2]
    b <- a + s$x[3]
    c <- (s$x[2] + s$x[3])*s$x[1] + s$x[2]*s$x[3]
            # coefficients for quadratic equation for minimum   
    A <- 3*s$q[4]                      # coeff of x^2 
    B <- (s$q[3] - s$q[4]*b)*2         # coeff of x
    C <- s$q[2] - a*s$q[3] + c*s$q[4]  # absolute term
  } else { #  Thiele
    z0 <- s$q[bm];  z1 <- z2 <- n1 <- n2 <- 0;    n0 <- 1;
    for (ii in ((bm-1):1)) {
      sq <- s$q[ii];  sx <- s$x[ii]
      if (abs(z0) < sqrt(sqrt(.Machine$double.xmax))) {
        Z0 <- sq*z0 - sx*n0
        Z1 <- sq*z1 - sx*n1 + n0
        Z2 <- sq*z2 - sx*n2 + n1
        # not needed    z3 <- n2  
        n0  <- z0;  n1 <- z1;  n2 <- z2
        z0  <- Z0;  z1 <- Z1;  z2 <- Z2
      } else {  # "restart"
        z0 <- sq;  z1 <- z2 <- n1 <- n2 <- 0;    n0 <- 1;
      }
    }
 # coefficients for quadratic equation for minimum   
    A <- z2*n1 - z1*n2      # coeff of x^2
    B <- 2*(z2*n0 - z0*n2)  # coeff of x
    C <- z1*n0 - z0*n1      # absolute term
  }
  xmin <- sort(solveQeq( A, B, C ), na.last=TRUE)
  if (Re(xmin[1]) <= min(x) || Re(xmin[1] >= max(x)) ) {
    xa <- xmin[1];  xmin[1] <- xmin[2];  xmin[2] <- xa 
  }
  if (!doPoly && abs(n2)/16.0+1.0 == 1.0) { # high chance of unreachable points
    if ( abs(n1)/16+1.0 == 1.0) { # constant denominator
      xmin <- sort(solveQeq( z2, z1, z0 ), na.last=TRUE)
    } else {
      xn <- sort(solveQeq( 0.0, n1, n0 ), na.last=TRUE)[1]
      if (any(abs(xmin-xn)/16+1.0 == 1.0)) { # common zeros
        xmin <- NA           # of numerator and denominator
      }
    }
  }
  return( xmin )
}  ## minInterp

quadmin <- function(x, y) { # Newton interpolation with 3 values
  ## around minimum/maximum value of y
  if (missing(y)) {y <- x[,2]; x <- x[,1]}
  m <- which.min(y)
  if (length(x) != length(y) | length(y) < 3 |
         m == 1 | m == length(y)) return(NA)
  d1 <- (y[m]-y[m-1])/(x[m]-x[m-1])
  d2 <- (y[m+1]-y[m])/(x[m+1]-x[m])
  e1 <- (d2-d1)/(x[m+1]-x[m-1])
  return((x[m-1]+x[m])/2 - d1/(2*e1))
}

