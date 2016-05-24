#' sparsevar: A package to estimate multivariate time series models (such as VAR and 
#' VECM), under the sparsity hypothesis. 
#' 
#' It performs the estimation of the matrices of the models using penalized 
#' least squares methods such as LASSO, SCAD and MCP. 
#'
#' @section svar functions:
#' \code{estimateVAR}, \code{estimateVECM}, \code{simulateVAR}, \code{createSparseMatrix},
#' \code{mcSimulations}, \code{plotMatrix}, \code{plotVAR}, \code{plotComparisonVAR},
#' \code{l2norm}, \code{l1norm}, \code{lInftyNorm}, \code{maxNorm}, \code{frobNorm},
#' \code{spectralRadius}, \code{spectralNorm}
#' 
#' @docType package
#' @name svar 
#'
NULL