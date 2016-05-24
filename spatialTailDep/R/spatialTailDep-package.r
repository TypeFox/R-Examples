#' spatialTailDep
#'
#' The package \code{spatialTailDep} provides functions implementing the pairwise M-estimator of 
#' parametric spatial tail dependence models for distributions attracted to a max-stable law. 
#' This is a rank-based estimator, constructed as the minimizer of the distance between a vector of 
#' integrals of parametric pairwise tail dependence functions and the vector of their empirical counterparts.
#' It is especially suited for high-dimensional data since it relies on bivariate margins only and, as a consequence of the rank-based
#' approach, the univariate marginal distributions need not be estimated. 
#' For a complete description of the pairwise M-estimator, see Einmahl et al. (2014).
#' 
#' Currently, this package allows for estimation of the Brown-Resnick process and the Gaussian extreme 
#' value process (usually known as the Smith model) in two-dimensional space. The main function of this package is 
#' \code{\link{Mestimator}}, but several other functions are exported as well: \code{\link{tailInt}}
#' returns the integral of a parametric bivariate stable tail dependence function over the unit square, \code{\link{tailIntEmp}}
#' returns the integral of the bivariate empirical stable tail dependence function over the unit square, and 
#' \code{\link{AsymVar}} returns the asymptotic variance matrix for a list of pairs of locations defined 
#' by the user. 
#' 
#' The function \code{\link{Mestimator}} combines these functions: first it computes a
#' pilot estimator based on the Euclidian distance between the integrals of the parametric and the empirical
#' stable tail dependence functions, then it calculates a weight matrix, defined as the inverse of the asymptotic
#' variance matrix in the pilot estimator, and finally it returns the estimator obtained by 
#' replacing the Euclidian distance by a quadratic form based on the weight matrix. More details
#' about this procedure can be found in Einmahl et al (2014).
#' 
#' The package exports two auxiliary functions as well: the function \code{\link{selectPairIndices}} returns
#' a list of pair indices of locactions, based on either a maximum-distance criterion or on the maximal 
#' number of pairs that the user wants to include. Next, the function \code{\link{pairCoordinates}} can be 
#' used to select the pairs with these indices from a list of location coordinates.  
#' 
#' @name spatialTailDep
#' @docType package
#' @references Einmahl, J.H.J., Kiriliouk, A., Krajina, A. and Segers, J., "An M-estimator of spatial tail dependence". See \url{http://arxiv.org/abs/1403.1975}. 
#' @examples
#' ## get a list of all help files of user-visible functions in the package
#' help(package = spatialTailDep)
NULL
  