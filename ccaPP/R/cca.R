# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' (Robust) CCA via alternating series of grid searches
#' 
#' Perform canoncial correlation analysis via projection pursuit based on 
#' alternating series of grid searches in two-dimensional subspaces of each 
#' data set, with a focus on robust and nonparametric methods.
#' 
#' The algorithm is based on alternating series of grid searches in 
#' two-dimensional subspaces of each data set.  In each grid search, 
#' \code{nGrid} grid points on the unit circle in the corresponding plane are 
#' obtained, and the directions from the center to each of the grid points are 
#' examined.  In the first iteration, equispaced grid points in the interval 
#' \eqn{[-\pi/2, \pi/2)}{[-pi/2, pi/2)} are used.  In each subsequent 
#' iteration, the angles are halved such that the interval 
#' \eqn{[-\pi/4, \pi/4)}{[-pi/4, pi/4)} is used in the second iteration and so 
#' on.  If only one data set is multivariate, the algorithm simplifies 
#' to iterative grid searches in two-dimensional subspaces of the corresponding 
#' data set.
#' 
#' In the basic algorithm, the order of the variables in a series of grid 
#' searches for each of the data sets is determined by the average absolute 
#' correlations with the variables of the respective other data set.  Since 
#' this requires to compute the full \eqn{(p \times q)}{(p x q)} matrix of 
#' absolute correlations, where \eqn{p} denotes the number of variables of 
#' \code{x} and \eqn{q} the number of variables of \code{y}, a faster 
#' modification is available as well.  In this modification, the average 
#' absolute correlations are computed over only a subset of the variables of 
#' the respective other data set.  It is thereby possible to use randomly 
#' selected subsets of variables, or to specify the subsets of variables 
#' directly.
#' 
#' Note that also the data sets are ordered according to the maximum average 
#' absolute correlation with the respective other data set to ensure symmetry 
#' of the algorithm.
#' 
#' For higher order canonical correlations, the data are first transformed into 
#' suitable subspaces.  Then the alternate grid algorithm is applied to the 
#' reduced data and the results are back-transformed to the original space.
#' 
#' @aliases print.cca
#' 
#' @param x,y  each can be a numeric vector, matrix or data frame.
#' @param k  an integer giving the number of canonical variables to compute.
#' @param method  a character string specifying the correlation functional to 
#' maximize.  Possible values are \code{"spearman"} for the Spearman 
#' correlation, \code{"kendall"} for the Kendall correlation, \code{"quadrant"} 
#' for the quadrant correlation, \code{"M"} for the correlation based on a 
#' bivariate M-estimator of location and scatter with a Huber loss function, or 
#' \code{"pearson"} for the classical Pearson correlation (see 
#' \code{\link{corFunctions}}).
#' @param control  a list of additional arguments to be passed to the specified 
#' correlation functional.  If supplied, this takes precedence over additional 
#' arguments supplied via the \code{\dots} argument.
#' @param nIterations,maxiter  an integer giving the maximum number of 
#' iterations.
#' @param nAlternate,maxalter  an integer giving the maximum number of 
#' alternate series of grid searches in each iteration.
#' @param nGrid,splitcircle  an integer giving the number of equally spaced 
#' grid points on the unit circle to use in each grid search.
#' @param select  optional; either an integer vector of length two or a list 
#' containing two index vectors.  In the first case, the first integer gives 
#' the number of variables of \code{x} to be randomly selected for determining 
#' the order of the variables of \code{y} in the corresponding series of grid 
#' searches, and vice versa for the second integer.  In the latter case, the 
#' first list element gives the indices of the variables of \code{x} to be used 
#' for determining the order of the variables of \code{y}, and vice versa for 
#' the second integer (see \dQuote{Details}).
#' @param tol,zero.tol  a small positive numeric value to be used for 
#' determining convergence.
#' @param standardize  a logical indicating whether the data should be 
#' (robustly) standardized.
#' @param fallback  logical indicating whether a fallback mode for robust 
#' standardization should be used.  If a correlation functional other than the 
#' Pearson correlation is maximized, the first attempt for standardizing the 
#' data is via median and MAD.  In the fallback mode, variables whose MADs are 
#' zero (e.g., dummy variables) are standardized via mean and standard 
#' deviation.  Note that if the Pearson correlation is maximized, 
#' standardization is always done via mean and standard deviation.
#' @param seed  optional initial seed for the random number generator (see 
#' \code{\link{.Random.seed}}).  This is only used if \code{select} specifies 
#' the numbers of variables of each data set to be randomly selected for 
#' determining the order of the variables of the respective other data set.
#' @param \dots  additional arguments to be passed to the specified correlation 
#' functional.  Currently, this is only relevant for the M-estimator.  For 
#' Spearman, Kendall and quadrant correlation, consistency at the normal model 
#' is always forced.
#' 
#' @returnClass cca
#' @returnItem cor  a numeric vector giving the canonical correlation 
#' measures.
#' @returnItem A  a numeric matrix in which the columns contain the canonical 
#' vectors for \code{x}.
#' @returnItem B  a numeric matrix in which the columns contain the canonical 
#' vectors for \code{y}.
#' @returnItem centerX  a numeric vector giving the center estimates used in 
#' standardization of \code{x}.
#' @returnItem centerY  a numeric vector giving the center estimates used in 
#' standardization of \code{y}.
#' @returnItem scaleX  a numeric vector giving the scale estimates used in 
#' standardization of \code{x}.
#' @returnItem scaleY  a numeric vector giving the scale estimates used in 
#' standardization of \code{y}.
#' @returnItem call  the matched function call.
#' 
#' @note \code{CCAgrid} is a simple wrapper function for \code{ccaGrid} for 
#' more compatibility with package \pkg{pcaPP} concerning function and argument 
#' names.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{ccaProj}}, \code{\link{maxCorGrid}}, 
#' \code{\link{corFunctions}}
#' 
#' @examples 
#' data("diabetes")
#' x <- diabetes$x
#' y <- diabetes$y
#' 
#' ## Spearman correlation
#' ccaGrid(x, y, method = "spearman")
#' 
#' ## Pearson correlation
#' ccaGrid(x, y, method = "pearson")
#' 
#' @keywords multivariate robust
#' 
#' @importFrom Rcpp evalCpp
#' @import robustbase
#' @useDynLib ccaPP
#' @export

ccaGrid <- function(x, y, k = 1, 
                    method = c("spearman", "kendall", "quadrant", "M", "pearson"), 
                    control = list(...), nIterations = 10, nAlternate = 10, 
                    nGrid = 25, select = NULL, tol = 1e-06, standardize = TRUE, 
                    fallback = FALSE, seed = NULL, ...) {
  ## initializations
  matchedCall <- match.call()
  ## define list of control arguments for algorithm
  nIterations <- as.integer(nIterations)
  nAlternate <- as.integer(nAlternate)
  nGrid <- as.integer(nGrid)
  tol <- as.numeric(tol)
  ppControl <- list(nIterations=nIterations, nAlternate=nAlternate, 
                    nGrid=nGrid, select=select, tol=tol)
  ## call workhorse function
  cca <- ccaPP(x, y, k, method=method, corControl=control, algorithm="grid", 
               ppControl=ppControl, standardize=standardize, fallback=fallback, 
               seed=seed)
  cca$call <- matchedCall
  cca
}

## wrapper function for more compatibility with package pcaPP
#' @rdname ccaGrid
#' @export

CCAgrid <- function(x, y, k = 1, 
                    method = c("spearman", "kendall", "quadrant", "M", "pearson"), 
                    maxiter = 10, maxalter = 10, splitcircle = 25, select=NULL, 
                    zero.tol = 1e-06, standardize = TRUE, fallback = FALSE, 
                    seed = NULL, ...) {
  ## initializations
  matchedCall <- match.call()
  ## call ccaGrid()
  cca <- ccaGrid(x, y, k=k, method=method, nIterations=maxiter, 
                 nAlternate=maxalter, nGrid=splitcircle, select=select, 
                 tol=zero.tol, standardize=standardize, fallback=fallback, 
                 seed=seed, ...)
  cca$call <- matchedCall
  cca
}


#' (Robust) CCA via projections through the data points
#' 
#' Perform canoncial correlation analysis via projection pursuit based on 
#' projections through the data points, with a focus on robust and 
#' nonparametric methods.
#' 
#' First the candidate projection directions are defined for each data set 
#' from the respective center through each data point.  Then the algorithm 
#' scans all \eqn{n^2} possible combinations for the maximum correlation, 
#' where \eqn{n} is the number of observations.
#' 
#' For higher order canonical correlations, the data are first transformed into 
#' suitable subspaces.  Then the alternate grid algorithm is applied to the 
#' reduced data and the results are back-transformed to the original space.
#' 
#' @param x,y  each can be a numeric vector, matrix or data frame.
#' @param k  an integer giving the number of canonical variables to compute.
#' @param method  a character string specifying the correlation functional to 
#' maximize.  Possible values are \code{"spearman"} for the Spearman 
#' correlation, \code{"kendall"} for the Kendall correlation, \code{"quadrant"} 
#' for the quadrant correlation, \code{"M"} for the correlation based on a 
#' bivariate M-estimator of location and scatter with a Huber loss function, or 
#' \code{"pearson"} for the classical Pearson correlation (see 
#' \code{\link{corFunctions}}).
#' @param control  a list of additional arguments to be passed to the specified 
#' correlation functional.  If supplied, this takes precedence over additional 
#' arguments supplied via the \code{\dots} argument.
#' @param standardize  a logical indicating whether the data should be 
#' (robustly) standardized.
#' @param useL1Median  a logical indicating whether the \eqn{L_{1}}{L1} medians 
#' should be used as the centers of the data sets in standardization (defaults 
#' to \code{TRUE}).  If \code{FALSE}, the columnwise centers are used instead 
#' (columnwise means if \code{method} is \code{"pearson"} and columnwise 
#' medians otherwise).
#' @param fallback  logical indicating whether a fallback mode for robust 
#' standardization should be used.  If a correlation functional other than the 
#' Pearson correlation is maximized, the first attempt for standardizing the 
#' data is via median and MAD.  In the fallback mode, variables whose MADs are 
#' zero (e.g., dummy variables) are standardized via mean and standard 
#' deviation.  Note that if the Pearson correlation is maximized, 
#' standardization is always done via mean and standard deviation.
#' @param \dots  additional arguments to be passed to the specified correlation 
#' functional.  Currently, this is only relevant for the M-estimator.  For 
#' Spearman, Kendall and quadrant correlation, consistency at the normal model 
#' is always forced.
#' 
#' @returnClass cca
#' @returnItem cor  a numeric vector giving the canonical correlation 
#' measures.
#' @returnItem A  a numeric matrix in which the columns contain the canonical 
#' vectors for \code{x}.
#' @returnItem B  a numeric matrix in which the columns contain the canonical 
#' vectors for \code{y}.
#' @returnItem centerX  a numeric vector giving the center estimates used in 
#' standardization of \code{x}.
#' @returnItem centerY  a numeric vector giving the center estimates used in 
#' standardization of \code{y}.
#' @returnItem scaleX  a numeric vector giving the scale estimates used in 
#' standardization of \code{x}.
#' @returnItem scaleY  a numeric vector giving the scale estimates used in 
#' standardization of \code{y}.
#' @returnItem call  the matched function call.
#' 
#' @note \code{CCAproj} is a simple wrapper function for \code{ccaProj} for 
#' more compatibility with package \pkg{pcaPP} concerning function names.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{ccaGrid}}, \code{\link{maxCorProj}}, 
#' \code{\link{corFunctions}}
#' 
#' @examples 
#' data("diabetes")
#' x <- diabetes$x
#' y <- diabetes$y
#' 
#' ## Spearman correlation
#' ccaProj(x, y, method = "spearman")
#' 
#' ## Pearson correlation
#' ccaProj(x, y, method = "pearson")
#' 
#' @keywords multivariate robust
#' 
#' @importFrom Rcpp evalCpp
#' @import pcaPP
#' @import robustbase
#' @useDynLib ccaPP
#' @export

ccaProj <- function(x, y, k = 1, 
                    method = c("spearman", "kendall", "quadrant", "M", "pearson"), 
                    control = list(...), standardize = TRUE, useL1Median = TRUE, 
                    fallback = FALSE, ...) {
  ## initializations
  matchedCall <- match.call()
  ## define list of control arguments for algorithm
  ppControl <- list(useL1Median=isTRUE(useL1Median))
  ## call workhorse function
  cca <- ccaPP(x, y, k, method=method, corControl=control, algorithm="proj", 
               ppControl=ppControl, standardize=standardize, fallback=fallback)
  cca$call <- matchedCall
  cca
}

## wrapper function for more compatibility with package pcaPP
#' @rdname ccaProj
#' @export

CCAproj <- function(x, y, k = 1, 
                    method = c("spearman", "kendall", "quadrant", "M", "pearson"), 
                    standardize = TRUE, useL1Median = TRUE, fallback = FALSE, 
                    ...) {
  ## initializations
  matchedCall <- match.call()
  ## call ccaProj()
  cca <- ccaProj(x, y, k=k, method=method, standardize=standardize, 
                 useL1Median=useL1Median, fallback=fallback, ...)
  cca$call <- matchedCall
  cca
}


## workhorse function
ccaPP <- function(x, y, k = 1, 
                  method = c("spearman", "kendall", "quadrant", "M", "pearson"), 
                  corControl, forceConsistency = TRUE, 
                  algorithm = c("grid", "proj"), ppControl, standardize = TRUE, 
                  fallback = FALSE, seed = NULL) {
  ## initializations
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  if(nrow(y) != n) {
    stop("'x' and 'y' must have the same number of observations")
  }
  p <- ncol(x)
  q <- ncol(y)
  # check number of canonical variables to compute
  k <- rep(as.integer(k), length.out=1)
  if(is.na(k) || k < 0) k <- formals()$k
  k <- min(k, p, q)
  ## prepare the data and call C++ function
  if(n == 0 || p == 0 || q == 0 || k == 0) {
    # zero dimension
    A <- B <- matrix(numeric(), 0, 0)
    cca <- list(cor=NA, A=A, B=B)
  } else {
    # check high-dimensional data
    if(k > 1 && (n <= p+1 || n <= q+1)) {
      k <- 1
      warning("higher-order canonical correlations not yet implemented", 
              "for high-dimensional data")
    }
    # check method and get list of control arguments
    method <- match.arg(method)
    corControl <- getCorControl(method, corControl, forceConsistency)
    # additional checks for grid search algorithm
    if(algorithm == "grid") {
      # check subset of variables to be used for determining the order of 
      # the variables from the respective other data set
      select <- ppControl$select
      ppControl$select <- NULL
      if(!is.null(select)) {
        if(is.list(select)) {
          # make sure select is a list with two index vectors and 
          # drop invalid indices from each vector
          select <- rep(select, length.out=2)
          select <- mapply(function(indices, max) {
            indices <- as.integer(indices)
            indices[which(indices > 0 & indices <= max)] - 1
          }, select, c(p, q))
          valid <- sapply(select, length) > 0
          # add the two index vectors to control object
          if(all(valid)) {
            ppControl$selectX <- select[[1]]
            ppControl$selectY <- select[[2]]
          } else select <- NULL
        } else {
          # check number of indices to sample
          select <- rep(as.integer(select), length.out=2)
          valid <- !is.na(select) & select > 0 & select < c(p, q)
          if(all(valid)) {
            # generate index vectors and add them to control object
            if(!is.null(seed)) set.seed(seed)
            ppControl$selectX <- sample.int(p, select[1]) - 1
            ppControl$selectY <- sample.int(q, select[2]) - 1
          } else select <- NULL
        }
      }
      if(is.null(select)) {
        ppControl$selectX <- ppControl$selectY <- integer()
      }
    }
    # call C++ function
    cca <- .Call("R_ccaPP", R_x=x, R_y=y, R_k=k, R_method=method, 
                 R_corControl=corControl, R_algorithm=algorithm, 
                 R_ppControl=ppControl, R_standardize=isTRUE(standardize), 
                 R_fallback=isTRUE(fallback), PACKAGE="ccaPP")
    drop <- c("cor", "centerX", "centerY", "scaleX", "scaleY")
    cca[drop] <- lapply(cca[drop], drop)
  }
  ## assign class and return results
  class(cca) <- "cca"
  cca
}
