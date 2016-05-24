##' Function to convert a correlation matrix to a covariance matrix.
##'
##' The correlation matrix to convert can be either symmetric or triangular. The covariance matrix returned is always a symmetric matrix.
##' @title Correlation Matrix to Covariance Matrix Conversion
##' @param cor.mat the correlation matrix to be converted
##' @param sd a vector that contains the standard deviations of the variables in the correlation matrix
##' @param discrepancy a neighborhood of 1, such that numbers on the main diagonal of the correlation matrix will be considered as equal to 1 if they fall in this neighborhood
##' @export
##' @note The correlation matrix input should be a square matrix, and the length of sd should be equal to the number of variables in the correlation matrix (i.e., the number of rows/columns). Sometimes the correlation matrix input may not have exactly 1's on the main diagonal, due to, eg, rounding; discrepancy specifies the allowable discrepancy so that the function still considers the input as a correlation matrix and can proceed (but the function does not change the numbers on the main diagonal). 
##' @author Ken Kelley (University of Notre Dame; (\email{KKelley@@ND.Edu}) and Keke Lai (the \code{MBESS} package), with modifications by Dustin Fife \email{fife.dustin@@gmail.com}.
cor2cov = function (cor.mat, sd, discrepancy = 0.00001) 
{
    if (dim(cor.mat)[1] != dim(cor.mat)[2]) 
        stop("'cor.mat' should be a square matrix")
    n <- sqrt(length(cor.mat))
    if (n != length(sd)) 
        stop("The length of 'sd' should be the same as the number of rows of 'cor.mat'")
    if (length(sd[sd > 0]) != n) 
        stop("The elements in 'sd' shuold all be positive")
    if (isSymmetric(cor.mat)) 
        IS.symmetric <- TRUE
    else IS.symmetric <- FALSE
    p <- dim(cor.mat)[1]
    q <- p * (p - 1)/2
    if (isTRUE(all.equal(cor.mat[lower.tri(cor.mat)], rep(0, 
        q))) || isTRUE(all.equal(cor.mat[upper.tri(cor.mat)], 
        rep(0, q)))) 
        IS.triangular <- TRUE
    else IS.triangular <- FALSE
    if (!IS.symmetric & !IS.triangular) 
        stop("The object 'cor.mat' should be either a symmetric or a triangular matrix")
    cov.mat <- diag(sd) %*% cor.mat %*% diag(sd)
    colnames(cov.mat) <- rownames(cov.mat) <- colnames(cor.mat)
    return(cov.mat)
}
