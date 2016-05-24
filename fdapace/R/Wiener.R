#' Simulate standard Wiener processes (Brownian motions)
#' 
#' Simulate \code{n} standard Wiener processes on [0, 1], possibly
#' sparsifying the results.
#'
#' The algorithm is based on Karhunen-Loeve expansion.
#' 
#' @param n Sample size.
#' @param pts A vector of points in [0, 1] specifying the support of the processes.
#' @param sparsify A vector of integers. The number of observations per curve will be uniform distribution on sparsify.
#' @param K The number of components.
#'
#' @return If \code{sparsify} is not specified, a matrix with \code{n} rows corresponding to the samples are returned. Otherwise the sparsified sample is returned. 
#'
#' @seealso Sparsify
#' @export

Wiener <- function(n=1, pts=seq(0, 1, length=50), sparsify=NULL, K=50) {
    pts <- as.matrix(pts)
    if (dim(pts)[1] < dim(pts)[2])
        pts <- t(pts)
        
    basis <- sqrt(2) * sin( pts %*% matrix(1:K - 1/2, 1, K) * pi )
    samp <- t(basis %*% diag(1 / (1:K - 1/2) / pi) %*% matrix(rnorm(K * n), K, n))
    
    if (!is.null(sparsify)) {
        samp <- Sparsify(samp, pts, sparsify)
    }
    
    return(samp)

}

## sparsify samp
## samp: a matrix of samples, with rows containing the samples
## pts: a vector of grid points, should be from 0 to 1
## sparsity: a vector of integers. The number of observation will be uniform distribution on sparsify.
#Sparsify <- function(samp, pts, sparsity) {
#    if (length(sparsity) == 1)
#            sparsity <- c(sparsity, sparsity) # avoid scaler case
#    
#    indEach <- lapply(1:nrow(samp), function(x) 
#        sort(sample(ncol(samp), sample(sparsity, 1))))
#    tList <- lapply(indEach, function(x) pts[x])
#    yList <- lapply(1:length(indEach), function(x) {
#        ind <- indEach[[x]]
#        y <- samp[x, ind]
#        return(y)
#    })
#   
#    return(list(tList=tList, yList=yList))
#}

