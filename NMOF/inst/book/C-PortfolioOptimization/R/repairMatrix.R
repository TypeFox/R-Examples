repairMatrix <- function(C) {
    # compute eigenvectors/-values
    E <- eigen(C, symmetric = TRUE)   
    V <- E$vectors
    D <- E$values
    
    # replace negative eigenvalues by zero
    D <- pmax(D, 0)
    
    # reconstruct correlation matrix
    BB <- V %*% diag(D) %*% t(V)
    
    # rescale correlation matrix
    T  <- 1/sqrt(diag(BB))
    TT <- outer(T,T)
    C  <- BB * TT
    C
}