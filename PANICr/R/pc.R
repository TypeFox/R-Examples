#'@title Principle Component Analysis for PANIC (2004)
#'
#'@description This function performs the principle component analysis in order to determine
#' the Common and Idiosyncratic Components of the factor model.
#'
#'@usage pc(y,nfac)
#'
#'
#'@param y A NxT matrix containing the data
#'
#'@param nfac An integer specifying the maximum number of factors allowed
#' while estimating the factor model.
#'
#'@return ehat The Idiosyncratic component of the factor model
#'
#'@return fhat The approximate factors of the approximate factor model
#'
#'@return lambda The factor loadings of the factor model
pc <- function(y, nfac) {
    
    bigt <- dim(y)[1]
    
    bign <- dim(y)[2]
    
    eig <- svd(crossprod(y))
    
    Fhat0 <- eig$u
    
    eigval <- as.matrix(eig$d)
    
    Fhat1 <- eig$v
    
    lambda <- Fhat0[, 1:nfac] * sqrt(bign)
    
    fhat <- y %*% (lambda/bign)
    
    ehat <- y - tcrossprod(fhat, lambda)
    
    return(list(ehat = ehat, fhat = fhat, lambda = lambda))
} 
