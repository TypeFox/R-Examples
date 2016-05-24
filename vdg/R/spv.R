#' Calculate the Scaled Prediction Variance (or SPV)
#' 
#' Calculates the SPV for a sample of points in a design region of specified type. Sampling is done
#' by calling \code{\link{sampler}}.
#' 
#' @param n number of samples to take
#' @param design a design or list of designs. Each design must be either a matrix or a data.frame or coercible to a data.frame.
#' @param type type of sampling passed to \code{\link{sampler}}
#' @param formula either a single model formula of a list of formulae
#' @param at only used when type is \code{'spherical'} or \code{'cuboidal'}
#' @param keepfun optional; function operating on the columns of a matrix with the same number of columns as design which return a logical value for 
#' including a specific point in the sample or not. Useful for rejection sampling for nonstandard design regions.
#' @param sample optional; if not missing it should contain a matrix or data.frame containing points sampled over the required design region. If it is not 
#' missing, no further sampling will be done: the SPV is simply evaluated at these points.
#' @param unscaled logical indicating whether to use the unscaled prediction variance (UPV) instead of the scale prediction variance (SPV)
#' @param \dots additional arguments passed to \code{\link{sampler}}. This enables the used of 
#' user-specified sampling functions via the \code{custom.fun} argument to \code{\link{sampler}}.
#' @return Object of class 'spv', 'spvlist', 'spvforlist' or 'spvlistforlist', depending on whether single designs/formulas
#' are passed or lists of these. 
#' @author Pieter C. Schoonees
#' @seealso \code{\link{plot.spv}} for more examples
#' @keywords multivariate
#' @export
#' @import parallel
#' @examples
#' 
#' # Single design (class 'spv')
#' library(rsm)
#' bbd3 <- as.data.frame(bbd(3)[,3:5])
#' colnames(bbd3) <- paste0("x", 1:3)
#' quad.3f <- formula(~(x1 + x2 + x3)^2 - x1:x2:x3)
#' out <- spv(n = 1000, design = bbd3, type = "spherical", formula = quad.3f)
#' out
#' @rdname spv
#' @export
spv <- function(n, design, type = "spherical", formula, at = FALSE, keepfun, sample, unscaled = FALSE, ...){
  UseMethod("spv", design)
}