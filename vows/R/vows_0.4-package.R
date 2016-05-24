#' Voxelwise semiparametrics
#' 
#' This package efficiently performs inference on a large set of parametric or
#' semiparametric regressions that are "parallel" in the sense that they have a
#' common design matrix.  The functions are inspired by neuroimaging
#' applications, where the parallel models pertain to a grid of brain locations
#' known as voxels.
#' 
#' Functions ending in ".mp" ("massively parallel") are designed for responses
#' in the form of a (wide) matrix; functions ending in "4d" take
#' four-dimensional response data (e.g., a set of images) and convert it to
#' matrix form so that the corresponding ".mp" function can be applied.
#' Examples include \code{\link{lm.mp}} and \code{\link{lm4d}} for ordinary
#' linear models, \code{\link{rlrt.mp}} and \code{\link{rlrt4d}} for restricted
#' likelihood ratio tests (RLRTs) of a parametric null hypothesis vs. a smooth
#' alternative, and \code{\link{semipar.mp}} and \code{\link{semipar4d}} for
#' smoothing (see Reiss et al., 2014).  Functions for interactive visualization
#' (\code{\link{rlrtpanel}} and \code{\link{funkpanel}}) are also provided.
#' 
#' @name vows-package
#' @aliases vows-package vows
#' @docType package
#' @author Philip Reiss \email{phil.reiss@@nyumc.org}, Yin-Hsiu Chen
#' \email{enjoychen0701@@gmail.com}, Lei Huang \email{huangracer@@gmail.com},
#' Lan Huo, Ruixin Tan and Rong Jiao \email{jiaorong007@@gmail.com}
#' 
#' Maintainer: Philip Reiss \email{phil.reiss@@nyumc.org}
#' @references Reiss, P. T., Huang, L., Chen, Y.-H., Huo, L., Tarpey, T., and
#' Mennes, M. (2014). Massively parallel nonparametric regression, with an
#' application to developmental brain mapping. \emph{Journal of Computational
#' and Graphical Statistics}, \emph{Journal of Computational and Graphical
#' Statistics}, 23(1), 232--248.
#' @keywords package
NULL

#' @import fda gamm4 rpanel
NULL



