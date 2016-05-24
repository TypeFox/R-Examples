# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' (Robust) maximum correlation via alternating series of grid searches
#' 
#' Compute the maximum correlation between two data sets via projection pursuit 
#' based on alternating series of grid searches in two-dimensional subspaces of 
#' each data set, with a focus on robust and nonparametric methods.
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
#' @aliases print.maxCor
#' 
#' @param x,y  each can be a numeric vector, matrix or data frame.
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
#' @param nIterations  an integer giving the maximum number of iterations.
#' @param nAlternate  an integer giving the maximum number of alternate series 
#' of grid searches in each iteration.
#' @param nGrid  an integer giving the number of equally spaced grid points on 
#' the unit circle to use in each grid search.
#' @param select  optional; either an integer vector of length two or a list 
#' containing two index vectors.  In the first case, the first integer gives 
#' the number of variables of \code{x} to be randomly selected for determining 
#' the order of the variables of \code{y} in the corresponding series of grid 
#' searches, and vice versa for the second integer.  In the latter case, the 
#' first list element gives the indices of the variables of \code{x} to be used 
#' for determining the order of the variables of \code{y}, and vice versa for 
#' the second integer (see \dQuote{Details}).
#' @param tol  a small positive numeric value to be used for determining 
#' convergence.
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
#' functional.
#' 
#' @returnClass maxCor
#' @returnItem cor  a numeric giving the maximum correlation estimate.
#' @returnItem a  numeric; the weighting vector for \code{x}.
#' @returnItem b  numeric; the weighting vector for \code{y}.
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
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{maxCorProj}}, \code{\link{ccaGrid}}, 
#' \code{\link{corFunctions}}
#' 
#' @references 
#' A. Alfons, C. Croux and P. Filzmoser (2016) Robust maximum association 
#' between data sets: The \R Package \pkg{ccaPP}.  \emph{Austrian Journal of 
#' Statistics}, \bold{45}(1), 71--79.
#' 
#' A. Alfons, C. Croux and P. Filzmoser (2016) Robust maximum association 
#' estimators.  \emph{Journal of the American Statistical Association}.  DOI 
#' 10.1080/01621459.2016.1148609.  In press.
#' 
#' @examples 
#' data("diabetes")
#' x <- diabetes$x
#' y <- diabetes$y
#' 
#' ## Spearman correlation
#' maxCorGrid(x, y, method = "spearman")
#' maxCorGrid(x, y, method = "spearman", consistent = TRUE)
#' 
#' ## Pearson correlation
#' maxCorGrid(x, y, method = "pearson")
#' 
#' @keywords multivariate robust
#' 
#' @importFrom Rcpp evalCpp
#' @useDynLib ccaPP
#' @export

maxCorGrid <- function(x, y, 
                       method = c("spearman", "kendall", "quadrant", "M", "pearson"), 
                       control = list(...), nIterations = 10, 
                       nAlternate = 10, nGrid = 25, select = NULL, 
                       tol = 1e-06, standardize = TRUE, fallback = FALSE, 
                       seed = NULL, ...) {
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
  maxCor <- maxCorPP(x, y, method=method, corControl=control, 
                     algorithm="grid", ppControl=ppControl, 
                     standardize=standardize, fallback=fallback, 
                     seed=seed)
  maxCor$call <- matchedCall
  maxCor
}


#' (Robust) maximum correlation via projections through the data points
#' 
#' Compute the maximum correlation between two data sets via projection pursuit 
#' based on projections through the data points, with a focus on robust and 
#' nonparametric methods.
#' 
#' First the candidate projection directions are defined for each data set 
#' from the respective center through each data point.  Then the algorithm 
#' scans all \eqn{n^2} possible combinations for the maximum correlation, 
#' where \eqn{n} is the number of observations.
#' 
#' @param x,y  each can be a numeric vector, matrix or data frame.
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
#' functional.
#' 
#' @returnClass maxCor
#' @returnItem cor  a numeric giving the maximum correlation estimate.
#' @returnItem a  numeric; the weighting vector for \code{x}.
#' @returnItem b  numeric; the weighting vector for \code{y}.
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
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{maxCorGrid}}, \code{\link{ccaProj}}, 
#' \code{\link{corFunctions}}, 
#' 
#' @examples 
#' data("diabetes")
#' x <- diabetes$x
#' y <- diabetes$y
#' 
#' ## Spearman correlation
#' maxCorProj(x, y, method = "spearman")
#' maxCorProj(x, y, method = "spearman", consistent = TRUE)
#' 
#' ## Pearson correlation
#' maxCorProj(x, y, method = "pearson")
#' 
#' @keywords multivariate robust
#' 
#' @importFrom Rcpp evalCpp
#' @import pcaPP
#' @useDynLib ccaPP
#' @export

maxCorProj <- function(x, y, 
                       method = c("spearman", "kendall", "quadrant", "M", "pearson"), 
                       control = list(...), standardize = TRUE, 
                       useL1Median = TRUE, fallback = FALSE, ...) {
  ## initializations
  matchedCall <- match.call()
  ## define list of control arguments for algorithm
  ppControl <- list(useL1Median=isTRUE(useL1Median))
  ## call workhorse function
  maxCor <- maxCorPP(x, y, method=method, corControl=control, 
                     algorithm="proj", ppControl=ppControl, 
                     standardize=standardize, fallback=fallback)
  maxCor$call <- matchedCall
  maxCor
}


## workhorse function
maxCorPP <- function(x, y, ...) {
  ## call workhorse function for canonical correlation analysis
  maxCor <- ccaPP(x, y, forceConsistency=FALSE, ...)
  ## modify object and return results
  maxCor <- list(cor=maxCor$cor, a=drop(maxCor$A), b=drop(maxCor$B), 
                 centerX=maxCor$centerX, centerY=maxCor$centerY, 
                 scaleX=maxCor$scaleX, scaleY=maxCor$scaleY)
  class(maxCor) <- "maxCor"
  maxCor
}
