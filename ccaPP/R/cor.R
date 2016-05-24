# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Fast implementations of (robust) correlation estimators
#' 
#' Estimate the correlation of two vectors via fast C++ implementations, with a 
#' focus on robust and nonparametric methods.
#' 
#' \code{corPearson} estimates the classical Pearson correlation.  
#' \code{corSpearman}, \code{corKendall} and \code{corQuadrant} estimate the 
#' Spearman, Kendall and quadrant correlation, respectively, which are 
#' nonparametric correlation measures that are somewhat more robust.  
#' \code{corM} estimates the correlation based on a bivariate M-estimator of 
#' location and scatter with a Huber loss function, which is sufficiently 
#' robust in the bivariate case, but loses robustness with increasing dimension.
#' 
#' The nonparametric correlation measures do not estimate the same population 
#' quantities as the Pearson correlation, the latter of which is consistent at 
#' the bivariate normal model.  Let \eqn{\rho}{rho} denote the population 
#' correlation at the normal model.  Then the Spearman correlation estimates 
#' \eqn{(6/\pi) \arcsin(\rho/2)}{(6/pi) arcsin(rho/2)}, while the Kendall and 
#' quadrant correlation estimate 
#' \eqn{(2/\pi) \arcsin(\rho)}{(2/pi) arcsin(rho)}.  Consistent estimates are 
#' thus easily obtained by taking the corresponding inverse expressions.
#' 
#' The Huber M-estimator, on the other hand, is consistent at the bivariate 
#' normal model.
#' 
#' @rdname corFunctions
#' @name corFunctions
#' 
#' @param x,y  numeric vectors.
#' @param consistent  a logical indicating whether a consistent estimate at the 
#' bivariate normal distribution should be returned (defaults to \code{FALSE}).
#' @param prob  numeric; probability for the quantile of the 
#' \eqn{\chi^{2}}{chi-squared} distribution to be used for tuning the Huber 
#' loss function (defaults to 0.9).
#' @param initial  a character string specifying the starting values for the 
#' Huber M-estimator.  For \code{"quadrant"} (the default), \code{"spearman"} 
#' or \code{"kendall"}, the consistent version of the respecive correlation 
#' measure is used together with the medians and MAD's.  For \code{"pearson"}, 
#' the Pearson correlation is used together with the means and standard 
#' deviations.
#' @param tol  a small positive numeric value to be used for determining 
#' convergence.
#' 
#' @return The respective correlation estimate.
#' 
#' @note 
#' The Kendall correlation uses a naive \eqn{n^2} implementation if 
#' \eqn{n < 30} and a fast \eqn{O(n \log(n))}{O(n log(n))} implementation for 
#' larger values, where \eqn{n} denotes the number of observations.
#' 
#' Functionality for removing observations with missing values is currently not 
#' implemented.
#' 
#' @author Andreas Alfons, \eqn{O(n \log(n))}{O(n log(n))} implementation of 
#' the Kendall correlation by David Simcha
#' 
#' @seealso \code{\link{ccaGrid}}, \code{\link{ccaProj}}, 
#' \code{\link[stats]{cor}}
#' 
#' @examples 
#' ## generate data
#' library("mvtnorm")
#' set.seed(1234)  # for reproducibility
#' sigma <- matrix(c(1, 0.6, 0.6, 1), 2, 2)
#' xy <- rmvnorm(100, sigma=sigma)
#' x <- xy[, 1]
#' y <- xy[, 2]
#' 
#' ## compute correlations
#' 
#' # Pearson correlation
#' corPearson(x, y)
#' 
#' # Spearman correlation
#' corSpearman(x, y)
#' corSpearman(x, y, consistent=TRUE)
#' 
#' # Kendall correlation
#' corKendall(x, y)
#' corKendall(x, y, consistent=TRUE)
#' 
#' # quadrant correlation
#' corQuadrant(x, y)
#' corQuadrant(x, y, consistent=TRUE)
#' 
#' # Huber M-estimator
#' corM(x, y)
#' 
#' @keywords multivariate robust
#' 
#' @importFrom Rcpp evalCpp
#' @useDynLib ccaPP

NULL


## barebones version of the Pearson correlation
#' @rdname corFunctions
#' @export
corPearson <- function(x, y) {
    # initializations
    x <- as.numeric(x)
    y <- as.numeric(y)
    n <- length(x)
    if(length(y) != n) stop("'x' and 'y' must have the same length")
    if(n == 0) return(NA)  # zero length vectors
    # call C++ function
    .Call("R_corPearson", R_x=x, R_y=y, PACKAGE="ccaPP")
}

## barebones version of the Spearman correlation
#' @rdname corFunctions
#' @export
corSpearman <- function(x, y, consistent = FALSE) {
    # initializations
    x <- as.numeric(x)
    y <- as.numeric(y)
    n <- length(x)
    if(length(y) != n) stop("'x' and 'y' must have the same length")
    if(n == 0) return(NA)  # zero length vectors
    # call C++ function
    .Call("R_corSpearman", R_x=x, R_y=y, R_consistent=isTRUE(consistent), 
        PACKAGE="ccaPP")
}

## barebones version of the Kendall correlation
#' @rdname corFunctions
#' @export
corKendall <- function(x, y, consistent = FALSE) {
    # initializations
    x <- as.numeric(x)
    y <- as.numeric(y)
    n <- length(x)
    if(length(y) != n) stop("'x' and 'y' must have the same length")
    if(n == 0) return(NA)  # zero length vectors
    # call C++ function
    .Call("R_corKendall", R_x=x, R_y=y, R_consistent=isTRUE(consistent), 
        PACKAGE="ccaPP")
}

## barebones version of the quadrant correlation
#' @rdname corFunctions
#' @export
corQuadrant <- function(x, y, consistent = FALSE) {
    # initializations
    x <- as.numeric(x)
    y <- as.numeric(y)
    n <- length(x)
    if(length(y) != n) stop("'x' and 'y' must have the same length")
    if(n == 0) return(NA)  # zero length vectors
    # call C++ function
    .Call("R_corQuadrant", R_x=x, R_y=y, R_consistent=isTRUE(consistent), 
        PACKAGE="ccaPP")
}

## barebones version of the Huber-type M-estimator
#' @rdname corFunctions
#' @export
corM <- function(x, y, prob = 0.9, 
        initial = c("quadrant", "spearman", "kendall", "pearson"), 
        tol = 1e-06) {
    # initializations
    x <- as.numeric(x)
    y <- as.numeric(y)
    n <- length(x)
    if(length(y) != n) stop("'x' and 'y' must have the same length")
    if(n == 0) return(NA)  # zero length vectors
    prob <- as.numeric(prob)
    initial <- match.arg(initial)
    tol <- as.numeric(tol)
    # call C++ function
    .Call("R_corM", R_x=x, R_y=y, R_prob=prob, R_initial=initial, R_tol=tol, 
        PACKAGE="ccaPP")
}
