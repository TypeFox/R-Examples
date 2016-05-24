# Computation of MEWMA steady-state pdf (multivariate mean monitoring)
mewma.psi <- function(l, cE, p, type="cond", hs=0, r=20) {
  if ( l<=0 | l>1 )		stop("l has to be between 0 and 1")
  if ( cE<=0 )			stop("threshold c has to be positive")
  if ( p<1 )			stop("wrong dimension parameter")
  if ( hs<0 )			stop("wrong head start value")
  if ( r<4 )			stop("resolution too small")

  itype <- pmatch(tolower(type), c("cond", "cycl")) - 1
  if ( is.na(itype) )		stop("wrong type of steady-state density")
  
  zeug <- .C("mewma_psi", as.double(l), as.double(cE), as.integer(p), as.integer(itype), as.double(hs), as.integer(r),
                          ans=double(length=3*r+1), PACKAGE="spc")$ans
        
  zahl <- zeug[1]      
  PSI  <- zeug[1:r + 1]                          
  w    <- zeug[1:r + r+1]
  z    <- zeug[1:r + 2*r+1]
  
  l2 <- ( (1-l)/l )^2 
  fchi <- function(u, a) 2*a * dchisq( u^2/l^2, p, ncp=l2*a^2)/l^2
  
  if ( itype == 0 ) psi <- Vectorize(function(x) sum( w * PSI * fchi(sqrt(x), z))/zahl, "x")
  
  if ( itype == 1 ) {
    if ( hs < 1e-9 )  psi <- Vectorize(function(x) dchisq( x/l^2, p)/l^2 / zahl + sum( w * PSI * fchi(sqrt(x), z) ), "x")
    if ( hs >= 1e-9 ) psi <- Vectorize(function(x) fchi(sqrt(x), hs) / zahl + sum( w * PSI * fchi(sqrt(x), z) ), "x")
  }  
  
  psi
}