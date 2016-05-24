################################################################################
#' Estimate several chaotic invariants using linear regression
#' @description
#' Several chaotic invariants are estimated by using linear regression. This function
#' provides a common interface for the estimate of all these parameters (see \code{\link{corrDim}}, 
#' \code{\link{dfa}} and \code{\link{maxLyapunov}} for examples).
#' @param x An object containing all the information needed for the estimate.
#' @param regression.range Range of values on the x-axis on which the regression is performed.
#' @param do.plot Logical value. If TRUE (default value), a plot of the regression is shown.
#' @param ... Additional parameters.
#' @return An estimate of the proper chaotic invariant.
#' @references H. Kantz  and T. Schreiber: Nonlinear Time series Analysis (Cambridge university press)
#' @author Constantino A. Garcia
#' @export estimate
estimate = function(x, regression.range, do.plot,...){
  UseMethod("estimate")
}
  



#' Get the embedding dimensions used for compute a chaotic invariant. 
#' @param x An object containing all the information needed for the estimate.
#' @return A numeric vector with the embedding dimensions used for compute a
#'  chaotic invariant. 
#' @references H. Kantz  and T. Schreiber: Nonlinear Time series Analysis (Cambridge university press)
#' @author Constantino A. Garcia
#' @export embeddingDims
embeddingDims = function(x){
  UseMethod("embeddingDims")
}


embeddingDims.default = function(x){
  return(x$embedding.dims)
}


#' Get the radius of the neighborhoods  used for the  computations of 
#' a chaotic invariant. 
#' @param x An object containing all the information needed for the estimate of
#' the chaotic invariant.
#' @return A numeric vector with the radius of the neighborhoods  used for the  computations of 
#' a chaotic invariant. 
#' @references H. Kantz  and T. Schreiber: Nonlinear Time series Analysis (Cambridge university press)
#' @author Constantino A. Garcia
#' @export radius
radius = function(x){
  UseMethod("radius")
}

radius.default = function(x){
  return (x$radius)
}


#' Get the order of the nonlinear chaotic invariant. 
#' @param x An object containing all the information needed for the estimate of
#' the chaotic invariant.
#' @return A numeric vector with the radius of the neighborhoods  used for the  computations of 
#' a chaotic invariant. 
#' @seealso \code{\link{corrDim}}, \code{\link{sampleEntropy}} 
#' @references H. Kantz  and T. Schreiber: Nonlinear Time series Analysis (Cambridge university press)
#' @author Constantino A. Garcia
#' @export nlOrder
nlOrder = function(x){
  UseMethod("nlOrder")
}

#' Plot local scaling exponents
#' @description Plots the local scaling exponents of the correlation sum or
#'  the average Shannon  information (when computing information dimension). 
#' @param x An object containing all the information needed for the estimate of
#' the chaotic invariant.
#' @param ... Additional graphical parameters.
#' @references H. Kantz  and T. Schreiber: Nonlinear Time series Analysis (Cambridge university press)
#' @author Constantino A. Garcia
#' @export plotLocalScalingExp
plotLocalScalingExp = function(x,...){
  UseMethod("plotLocalScalingExp")
}
