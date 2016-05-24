grnmol <- function( f, V, s ){
#
#   R translation of Alan Genz's grnmol.m file for Grundmann-Moller integration
#   Translated by John Nolan, jpnolan@american.edu, 21 Aug 2014
#
#   [Q nv] = grnmol( f, V, s )
#     computes approximation to the integral of f over a simplex
#     with vertices as columns of V, an n x (n+1) matrix, using
#     order 1, 2, ..., s (degree 2s+1) Grundmann-Moler rules.
#   Output Q is a vector approximations of degree 1, 3, ... 2s+1. 
#   matlab example:  final two results should be 2/(11x10x9x8x7x6) ~ 6.012506e-06
#    n = 4; disp( grnmol(@(x)x(1)^2*x(n)^5,eye(n,n+1),4) )
#   R version of the above example:
#     n <- 4; V <- cbind( diag(rep(1.0,n)), rep(0.0,n) ); s <- 4
#     f <- function( x ) { x[1]^2*x[n]^5 }
#
#     Reference:
#       "Invariant Integration Formulas for the N-Simplex by
#        Combinatorial Methods", A. Grundmann and H. M. Moller, 
#        SIAM J Numer. Anal. 15(1978), pp. 282-290
#

# check input values
stopifnot( is.function(f), is.matrix(V), is.numeric(s), length(s)==1, s > 0, nrow(V)+1==ncol(V) )
s <- as.integer(s)

n <- nrow(V); np <- n+1
nd <- floor(n/2); d <- 0; Q <- rep(0.0,s+1); Qv <- Q;
Vol <- abs( det(V[,1:n]-V[,np] %*% matrix(1.0,nrow=1,ncol=n)) )/factorial(n)
repeat {
  m <- n + 2*d + 1; al <- rep(1.0,n); alz <- 2*d + 1; Qs <- 0; nv <- 0; 
  repeat {
    Qs <- Qs + f( V %*% c(alz,al)/m ); nv <- nv + 1;
    for (j in 1:n) {
      alz <- alz - 2; 
      if (alz > 0) { al[j] <- al[j] + 2; break }
      alz <- alz + al[j] + 1; al[j] <- 1;
    }
    if (alz == 2*d+1) { break }
  } # end of inner repeat loop
  d <- d + 1; Qv[d] <- Vol*Qs; Q[d] <- 0; p <- 2/prod( 2*((n+1):m) ); 
  for (i in 1:d) {
    Q[d] <- Q[d] + (m+2-2*i)^(2*d-1)*p*Qv[d+1-i]; p <- -p*(m+1-i)/i; 
  }
  if (d > s) { break }
} # end of outer repeat loop

return( list(Q=Q,nv=nv) ) }
# end grnmol

