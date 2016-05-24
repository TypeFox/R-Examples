# =========== Inner estimators ===========

#'@title PLS inner estimation 
#'
#'@description
#'
#'Calculates a set of inner weights. 
#
#'@details
#'In the centroid scheme, inner weights are set to the signs (1 or -1) of correlations between
#'composites that are connected in the model specified in \code{inner.mod} and zero otherwise.
#'
#'In the path scheme, inner weights are based on regression estimates of the relationships between
#'composites that are connected in the model specified in \code{inner.mod}, and correlations for
#'the inverse relationships. If a relationship is reciprocal, regression is used for both directions.
#'
#'In the factor scheme, inner weights are set to the correlations between
#'composites that are connected in the model specified in \code{inner.mod} and zero otherwise.
#'
#'In the identity scheme identity matrix is used as the inner weight matrix \code{E}.
#'
#'Centroid, innner, and path schemes fall back to to identity scheme for composites 
#'that are not connected to any other composites.
#'
#'For information about GSCA weights, see \link{GSCA}. 
#'
#'@inheritParams matrixpls-common
#'
#'@param inner.mod A square matrix specifying the relationships of the composites in the model.
#'
#'@param ignoreInnerModel Should the inner model be ignored and all correlations be used.
#'
#'@param ... Other arguments are ignored.
#'
#'@return A matrix of unscaled inner weights \code{E} with the same dimesions as \code{inner.mod}.
#'
#'
#'@references
#'Lohmöller J.-B. (1989) \emph{Latent variable path modeling with partial least squares.} Heidelberg: Physica-Verlag.
#'
#'@name innerEstimators
NULL

#'@title Generalized structured component analysis (GSCA) weights
#'
#'@description
#'
#'When used with \code{\link{weight.pls}}, \code{\link{inner.gsca}} and 
#'\code{\link{outer.gsca}} implement the generalized structured component analysis
#'indicator weighting system. Using \code{\link{weight.optim}} with the
#'\code{\link{optim.gsca}} optimization criterion provides an anternative
#'approach to calculate GSCA weights by direct numerical minimization of the
#'GSCA criterion function.
#'
#'@details
#'
#'The two step GSCA weight algorithm is designed to minimize.
#'
#'\code{SS(ZV-ZWA)}
#'
#'The parameter matrix \code{A} contains all model parameters including
#'\code{inner}, reflective \code{inner}, and \code{formative}. The weight
#'matrices \code{V} and \code{W}, which can contain duplicate elements,
#'contain the indicator weights.
#'
#'The first step of GSCA estimation method is calculation of regressions
#'coefficients \code{A} given weights \code{W} and \code{V}. The function
#'\code{\link{inner.gsca}} update the part of \code{A} corresponding to 
#'regressions between the composites, corresponding to \code{E} matrix in 
#'matrixpls. The regressions between composites and indicators are estimated
#'in \code{\link{outer.gsca}}.
#'
#'This algorithm for estimating the relationships between the composites
#'is therefore identical to the PLS path weighting scheme with
#'the exception that correlations are not used for inverse relationships and
#'there is no falling back to identity scheme for composites that are not
#'connected to other composites

#'The second step of GSCA is calculating a new set of weights conditional on
#'the regression coeffcients \code{A} to minimize the sum of error terms in
#'the regressions. This step is implemented in \code{\link{outer.gsca}} after
#'updating the regresions between indicators and composites.

#'The implementation of GSCA in \pkg{matrixpls} differs from the Hwang & Takane (2004)
#'version in that during the first step, only regressions between composites are
#'estimated. The regressions between composites and indicators are estimated by
#'the second stage 
#'the indicators and compositess. Since these covariances need to be calculated in the second step, it is more
#'efficient to not calculate them during the first step.
#'
#'A part of the code for \code{\link{outer.gsca}} is adopted from the \pkg{ASGCA} package, licenced
#'under GPL3.
#'
#'@author Mikko Rönkkö, Hela Romdhani, Stepan Grinek, Heungsun Hwang, Aurelie Labbe.
#'
#'@references
#'Hwang, H., & Takane, Y. (2004). Generalized structured component analysis. Psychometrika, 69(1), 81–99. doi:10.1007/BF02295841
#'
#'Hela Romdhani, Stepan Grinek, Heungsun Hwang and Aurelie Labbe. (2014). ASGSCA: Association Studies for multiple SNPs and multiple traits using
#' Generalized Structured Equation Models. R package version 1.4.0.
#' 
#'@example example/fragment-requireASGSCA.R
#'@example example/gsca-example.R
#'@example example/fragment-endBlock.R
#'
#'@name GSCA
#'
NULL


#'@describeIn innerEstimators inner estimation with centroid scheme.
#'@export

inner.centroid <- function(S, W, inner.mod, ignoreInnerModel = FALSE, ...){
  
  # Centroid is just the sign of factor weighting
  
  E <- sign(inner.factor(S, W, inner.mod, ignoreInnerModel, ...))
  
  return(E)
}

#'@describeIn innerEstimators inner estimation with path scheme.
#'@export

inner.path <- function(S, W, inner.mod, ...){
  
  # Calculate the composite covariance matrix
  C <- W %*% S %*% t(W)
  
  E <- estimator.ols(S, inner.mod,W, C = C)
  
  # Use correlations for inverse relationships for non-reciprocal paths
  inverseRelationships <- t(inner.mod) & ! inner.mod
  E[inverseRelationships] <- C[inverseRelationships]
  
  # If we have LVs that are not connected to any other LVs, use identity scheme as fallback
  diag(E)[rowSums(E) == 0] <- 1
  
  return(E)
}

#'@describeIn innerEstimators inner estimation with factor scheme.
#'@export

inner.factor <- function(S, W, inner.mod, ignoreInnerModel = FALSE, ...){
  
  # Calculate the composite covariance matrix
  
  C <- W %*% S %*% t(W)
  
  # Calculate unscaled inner weights using the factor weighting scheme
  if(ignoreInnerModel){
    E <- C
    diag(E) <- 0
  } 
  else E <- C * (inner.mod | t(inner.mod))
  
  # If we have LVs that are not connected to any other LVs, use identity scheme as fallback
  diag(E)[rowSums(E) == 0] <- 1
  
  return(E)
} 



#'@describeIn innerEstimators inner estimation with identity scheme.
#'@export

inner.identity <- function(S, W, inner.mod, ...){
  return(diag(nrow(inner.mod)))
}

#'@describeIn innerEstimators inner estimation with generalized structured component analysis.
#@describeIn GSCA inner estimation with generalized structured component analysis.
#'@export

inner.gsca <- function(S, W, inner.mod, ...){
  E <- estimator.ols(S, inner.mod,W)
  return(E)
}

