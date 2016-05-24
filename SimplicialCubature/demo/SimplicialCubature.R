# misc. examples for integrating over simplices
################################################################
# integrate over a m dim. surface in n dim.
fp1 <- function( x ) { x[1]^2 * x[2]^3 + 3}
for (n in 2:5) {
  cat("-------------------------------------------------------\n")
  cat("n m  integPoly/adaptInteg = ratio \n")
  for (m in 1:n) {
    S1 <- CanonicalSimplex( n=n )[ ,1:(m+1)]
    p1 <- definePoly( c(1.0,3.0), matrix( c(2,3,rep(0,2*n-2)), nrow=2, byrow=TRUE) )
    a <- integrateSimplexPolynomial( p1, S1, method="GM" ) 
    b <- adaptIntegrateSimplex( fp1, S1, maxEvals=1000000 )$integral
    cat(n,m,"   ",a$integral,"/",b,"=",a$integral/b,"\n")
  }
}
######################################################################
# example from integrateSimplexPolynomial
# full dimensional simplex
for (n in 2:5) {
  S <- CanonicalSimplex( n ) # 4-dim. simplex
  p1 <- definePoly( 1.0, matrix( c(2,3,rep(0,n-2)), nrow=1) )
  if(n==2) printPoly(p1)
  a1 <- integrateSimplexPolynomial( p1, S, method="GM" )
  a2 <- integrateSimplexPolynomial( p1, S, method="LA" )
  fp1 <- function( x ) { x[1]^2 * x[2]^3} 
  a3 <- adaptIntegrateSimplex( fp1, S )
  cat("n=",n,"   ",a1$integral,a2$integral,a3$integral,"\n")
}

