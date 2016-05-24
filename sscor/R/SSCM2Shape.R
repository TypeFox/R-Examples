# Calculates theoretic Shape matrix from theoretic SSCM
# input: theoretic SSCM
# tol: absolute tolerance for the approximated eigenvalues
# itermax: maximal number of iterations of the approximation algorithm 
# output: shape matrix

SSCM2Shape <- function(V,itermax=100,tol=10^(-10)) {
eigens <- eigen(V,symmetric=TRUE)
eigenc <- evSSCM2evShape(eigens$values,itermax=itermax,tol=tol)
return(eigens$vectors%*%diag(eigenc)%*%t(eigens$vectors))
}
