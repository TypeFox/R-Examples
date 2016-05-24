#' Graphical model check
#'
#' A graphical model check is performed for the multidimensional polytomous
#' Rasch model or the continuous Rating Scale Model.
#'
#' The graphical model check plots the item parameter estimates of two
#' subsamples to check the homogeneity. This is according to the subsample
#' split in Andersen's Likelihood Ratio test. For conducting the graphical
#' model check of the MPRM, at first, a \code{\link{LRT}} has to be computed
#' and the resulting object is the input for the \code{gmc} function.
#'
#' For plotting a graphical model check for the CRSM, the model has to be
#' estimated with \code{\link{CRSM}} and subsequently the resulting object is
#' the input for the \code{gmc} function. For the CRSM a split criterion has to
#' be input as vector.
#'
#' @aliases gmc gmc.aLR gmc.CRSM
#' @param object Object of class \code{aLR} for graphical model check of the
#' MPRM or object of class \code{CRSM} for graphical model check of the CRSM
#' @param \dots \dots{}
#' @param splitcrit Vector or the character vector \code{"score"} to define the
#' split criterion. The default split criterion \code{"score"} splits the
#' sample according to the median of the raw score. Vector can be numeric,
#' factor or character. (see details)
#' @author Christine Hohensinn
#' @seealso \code{\link{LRT}} \code{\link{CRSM}}
#' @references Wright, B.D., and Stone, M.H. (1999). Measurement Essentials.
#' Wilmington: Wide Range Inc.
#'
#' @rdname gmc
#' @keywords graphical model check
#' @examples
#'
#' #estimate CRSM for the first three items
#' data(analog)
#' res_cr <- CRSM(extraversion, low=-10, high=10)
#'
#' #graphical model check for CRSM for the first three items with default split
#' #criterion score
#' gmc(res_cr)
#'
#'
#' @export gmc
gmc <-
function(object,...)UseMethod("gmc")
