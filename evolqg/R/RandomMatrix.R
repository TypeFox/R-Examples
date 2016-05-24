#' Random matrices for tests
#'
#' Provides random covariance/correlation matrices for quick tests.
#' Should not be used for statistics or hypothesis testing.
#' @param num.traits Number of traits in random matrix
#' @param num.matrices Number of matrices to be generated. If greater than 1, a list is returned.
#' @param min.var Lower value for random variance in covariance matrices
#' @param max.var Upper value for random variance in covariance matrices
#' @param variance Variance vector. If present will be used in all matrices
#' @param ke Parameter for correlation matrix generation. Involves check for positive definiteness
#' @export
#' @return Returns either a single matrix, or a list of matrices of equal dimension
#' @author Diogo Melo Edgar Zanella
#' @examples
#' #single 10x10 correlation matrix
#' RandomMatrix(10)
#'
#' #single 5x5 covariance matrix, variances between 3 and 4
#' RandomMatrix(5, 1, 3, 4)
#'
#' #two 3x3 covariance matrices, with shared variances
#' RandomMatrix(3, 2, variance= c(3, 4, 5))
#'
#' #large 10x10 matrix list, with wide range of variances
#' RandomMatrix(10, 100, 1, 300)
#' @keywords randommatrices
RandomMatrix <- function(num.traits, num.matrices = 1, min.var = 1, max.var = 1, variance = NULL, ke = 10^-3){
    if(num.matrices==1){
        if(is.null(variance)) variance <- runif(num.traits, min.var, max.var)
        rand.mat <- RandCorr(num.traits, ke) * sqrt(outer(variance, variance))
    } else{
        if(is.null(variance)) variance <- matrix(runif(num.matrices*num.traits, min.var, max.var), num.matrices, num.traits)
        else variance <-  matrix(rep(variance, each = num.matrices), num.matrices, num.traits)
        rand.mat <- lapply(as.list(1:num.matrices), function(x) RandCorr(num.traits, ke) * sqrt(outer(variance[x,], variance[x,])))
    }
    return(rand.mat)
}
