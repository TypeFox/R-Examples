# Computation of GRSR (Girshick, Rubin, Shiryaev, Roberts) alarm threshold for given ARL (mean monitoring)
xgrsr.crit <- function(k, L0, mu0=0, zr=0, hs=NULL, sided="one", MPT=FALSE, r=30) {
  if ( k<0 )	stop("k has to be non-negative")
  if ( L0<1 )	stop("L0 is too small")
  if ( !is.null(hs) ) {
    if ( hs>log(L0) )	stop("wrong headstart")
  } else {
    hs <- 2*L0
  }
  if ( r<4 )	stop("r is too small")
  
  g <- .C("xgrsr_crit",as.double(k),
          as.double(L0),as.double(zr),as.double(hs),as.double(mu0),as.integer(r),as.integer(MPT),
          ans=double(length=1),PACKAGE="spc")$ans 
          
  names(g) <- "g"
  return (g)
}

