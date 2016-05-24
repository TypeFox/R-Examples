#' Function to run basic IPFP (iterative proportional fitting procedure)
#'
#' Use IPFP starting from x0 to produce vector x s.t. Ax = y within tolerance.
#' Need to ensure that x0 >= 0.
#'
#' @param y numeric constraint vector (length nrow)
#' @param A constraint matrix (nrow x ncol)
#' @param x0 numeric initial vector (length ncol)
#' @param tol numeric tolerance for IPFP; defaults to \code{.Machine$double.eps}
#' @param maxit integer maximum number of iterations for IPFP; defaults to 1e3
#' @param verbose logical parameter to select verbose output from C function
#' @param full logical parameter to select full return (with diagnostic info)
#' @return if not full, vector of length ncol containing solution obtained by
#'      IPFP. If full, list containing solution (as x), number of iterations (as
#'      iter), and norm of Ax - y (as errNorm)
#' @keywords iteration array
#' @export
#' @useDynLib networkTomography
#' @examples
#' A <- buildStarMat(3)
#' x <- rgamma(ncol(A), 10, 1/100)
#' y <- A %*% x
#' x0 <- x * rgamma(length(x), 10, 10)
#' ans <- ipfp(y, A, x0, full=TRUE)
#' print(ans)
#' print(x)
ipfp <- function(y, A, x0,
    tol=.Machine$double.eps, maxit=1e3, verbose=FALSE, full=FALSE) {
    # Get active rows
    activeRows <- which(y > 0)
    
    # Zero inactive columns
    if ( any(y==0) ) {
        activeCols <- !pmin(1, colSums(A[y==0,,drop=FALSE]))
    } else {
        activeCols <- rep(TRUE, ncol(A))
    }
    x0[!activeCols] <- 0
    x0[activeCols] <- pmax(1, x0[activeCols])
    
    # Run IPF
    ans <- .Call("ipfp", y[activeRows], A[activeRows, activeCols, drop=FALSE],
            dim(A[activeRows, activeCols, drop=FALSE]), x0[activeCols],
            as.numeric(tol), as.integer(maxit), as.logical(verbose),
            PACKAGE='networkTomography')
    
    x0[activeCols] <- ans$x
    
    if (full)
        return( list(x=x0, iter=ans$iter, errNorm=ans$errNorm) )
    return(x0)
}
