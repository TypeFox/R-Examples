#' Regression discontinuity estimation package
#' 
#' \code{rdd} supports both sharp and fuzzy RDD
#' utilizing the \pkg{AER} package for 2SLS regression
#' under the fuzzy design. Local linear regressions are performed
#' to either side of the cutpoint using the Imbens-Kalyanamaran
#' optimal bandwidth calculation, \code{\link{IKbandwidth}}.
#' 
#' @seealso \code{\link{RDestimate}}, \code{\link{DCdensity}}, \code{\link{IKbandwidth}},
#' \code{\link{summary.RD}}\code{\link{plot.RD}}, \code{\link{kernelwts}}
#' @name rdd-package
#' @aliases rdd
#' @docType package
#' @title Regression Discontinuity Estimation Package
#' @author Drew Dimmery \email{drewd@@nyu.edu}
NULL