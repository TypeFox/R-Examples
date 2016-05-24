#' Whitening Data Using Singular Value Decomposition
#' 
#' The function whitens a data matrix using the singular value decomposition.
#' 
#' The function whitens the data so that the result has mean zero and identity
#' covariance matrix using the function \code{\link{svd}}. The data can have
#' here less observations than variables and svd will determine the rank of the
#' data automatically as the number of singular values larger than the largest
#' singular value times \code{tol}. If \code{tol=NULL} the rank is set to the
#' number of singular values, which is not advised when one or more singular
#' values are close to zero.
#' 
#' The output contains among others as attributes the singular values and the
#' matrix needed to backtransform the whitened data to its original space.
#' 
#' @param x A numeric data frame or data matrix with at least two rows.
#' @param tol Tolerance value to decide the rank of the data. See details for
#' further information. If set to \code{NULL} it will be ignored.
#' @return The whitened data.
#' @author Klaus Nordhausen
#' @seealso \code{\link{svd}}, \code{\link{cov}}, \code{\link{colMeans}}
#' @keywords multivariate
#' @examples
#' 
#' # more observations than variables
#' X <- matrix(rnorm(200),ncol=4)
#' A <- WhitenSVD(X)
#' round(colMeans(A),4)
#' round(cov(A),4)
#' # how to backtransform
#' summary(sweep(A %*% (attr(A,"backtransform")), 2, attr(A,"center"), "+") - X)
#' 
#' # fewer observations than variables
#' Y <- cbind(rexp(4),runif(4),rnorm(4), runif(4), rnorm(4), rt(4,4))
#' B <- WhitenSVD(Y)
#' round(colMeans(B),4)
#' round(cov(B),4)
#' # how to backtransform
#' summary(sweep(B %*% (attr(B,"backtransform")), 2, attr(B,"center"), "+") - Y)
#' 
#' 
#' @export WhitenSVD
WhitenSVD <- function(x, tol=1e-06)
    {

    Nm <- nrow(x)-1
    MEANS <- colMeans(x)
    x.c <- as.matrix(sweep(x,2,MEANS,"-"))

    SVD <- svd(x.c, nu = 0)
    SV <- SVD$d 
    if (!is.null(tol)) {
        rank <- sum(SVD$d > (SVD$d[1L] * tol))
        if (rank < ncol(x)) {
            SVD$v <- SVD$v[, 1L:rank, drop = FALSE]
            SVD$d <- SVD$d[1L:rank]
            }
        }
    SIGMAs <- SVD$d / sqrt(Nm)
    TRANS <- SVD$v %*% diag(1/SIGMAs)
    RES = x.c %*% TRANS


    attr(RES, "center") <- MEANS
    attr(RES, "transform") <- TRANS
    attr(RES, "backtransform") <- diag(SIGMAs) %*% t(SVD$v)
    attr(RES, "SingularValues") <- SV
    RES
    }
