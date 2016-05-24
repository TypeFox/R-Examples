## nmlz: Returns the normalization matrix "J" of columns of M. That is: M %*% J has unit length columns. Zero columns correspond to "1" (unchanged).
nmlz = function(M){
    J = diag(1/apply(M,2,l2))
    diag(J)[which(diag(J) == Inf)] = 1
    return(J)
}
