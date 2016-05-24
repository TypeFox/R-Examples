lV <- function( v ) as.vector( sqrt(v %*% v) )

scprod <- function(v, w) { # v . w
  return( sum(v*w) )
} # end scprod

vecprod <- function(v, w) { # v cross w
  return( c(v[2]*w[3]  - v[3]*w[2],
            v[3]*w[1]  - v[1]*w[3],
            v[1]*w[2]  - v[2]*w[1] ))
} # end vecprod

"%s%" <- function(v,w) scprod(v,w)
"%v%" <- function(v,w) vecprod(v,w)

angle <- function(v, w)  {
  cs <- v %*% w
  lp <- lV(v)*lV(w)
  pc <- cs/lp
  if (abs(pc - 1) >= .Machine$double.eps*128.0 ) {
    res <- acos( pc )
  } else {  # small angles use asin instead of acos
    vp <- lV(vecprod(v, w))
    res <- asin(vp/lp)
  }
  return( res )
}

toPol <- function( x, y=0 ) {
  z <- complex( real=x,imaginary=y )
  res <- c( Mod(z), Arg(z))
  return ( res )
} ## end  toPol

toRec <- function( r, phi=0 ) {
  z <- complex(modulus=r,argument=phi)
  return ( c( Re(z), Im(z) ) )
}  ## end  toRect

toSph <- function ( x,y,z )  # inverse of toXyz
{
  switch( nargs() ,
    "1" = {z <-  x[3]; y <-  x[2]; x <- x[1]},
    "2" = {DG <-  x[4]; z <-  x[3]; y <-  x[2]; x <- x[1]},
    "3" = {},
    "4" = {}
         )
  vphi <- toPol( x,y ) # x, y -> c( w, phi )
  R <- toPol( z, vphi[1] )  # ( w, z,  -> r, theta )
  res <- c(R[1], R[2], vphi[2]) # r, theta, phi
  return ( res )
} ## end  toSph

toXyz <-  function ( r, theta, phi )  # inverse of toSph
{
  switch( nargs() ,
    "1" = {phi <-  r[3]; theta <-  r[2]; r <- r[1]},
    "2" = {DG <-  r[4]; phi <-  r[3]; theta <-  r[2]; r <- r[1]},
    "3" = {},
    "4" = {}
         )
  vz <- toRec( r, theta ) # ->  z, w
  xy <- toRec( vz[2], phi )
  res <- c(xy, vz[1] ) # x, y, z
  return ( res )
} ## end  toXyz

rotZ <- function(x, y, phi) {
  rp <- toPol( x, y )
  return ( toRec(rp[1], rp[2] + phi) )
}  ## end rotZ

rotA <- function( phi, P=c(0,0,1) )  {
##  r <- rotC( P ) ## rotate P into z-axis
  r <- rotL( -acos(P[3]/sqrt(t(P)%*%P)) ,1,3) %*% rotL( atan2(P[2],P[1]) ,1,2)
  t(r) %*% rotL( phi, 1, 2 ) %*% r
}

rotV <- function(v, w=c(0,0,1)) {
  u <- vecprod(v, w)
  if ( lV(u) <= .Machine$double.eps*16.0) {
    res <- diag(3)  # v, w almost parallel
  } else {
    phi <- angle(v,w)
    res <- rotA(phi, u)
  }
  return( res )
}

rotL <- function(phi,k=1,m=2,N=3) {
  res <- diag(N)
  if (k != m) {
    ss  <- sin(phi)
    cc  <- cos(phi)
    res[k,k] <- res[m,m] <- cc
    res[k,m] <- ss
    res[m,k] <- -ss
  }
  return( res )
}

getAp <- function( M ) { # determine axis and angle from rotation matrix
  e <- eigen(M)
  return(list(A=Re(e$vectors[,1]),phi=acos(Re(sum(e$values)-1)/2)))
}

