#
# vim:set ff=unix expandtab ts=2 sw=2:
spectralNorm=function(#a utility function to compute the spectral norm of a matrix
### This function computes the spectral norm of a matrix which is defined as the inverse of the minimum of absolute values of the eigenvalues.
m ##<<a (quadratic) matrix 
){

1/min(abs((eigen(m,only.values=TRUE))$values))
### The spectral norm of the matrix, a real number.
}

