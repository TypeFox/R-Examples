#' Item Characteristic Curve
#'
#' The item characteristic curve is performed for the multidimensional polytomous
#' Rasch model or the continuous Rating Scale Model.
#'
#' The item characteristic curve (ICC) plots the response probability depending on person and item parameter.
#' For plotting the ICC, the object resulting from MPRM \code{\link{MPRM}} or CRSM \code{\link{CRSM}} or DRM \code{\link{DRM}} is the input for the \code{iccplot} function.
#' The default argument \code{items="all"} displays ICC curves for all items in the object. With a numeric vector \code{items}, a subset of
#' items can be selected for which ICC plots are displayed.
#'
#' @aliases iccplot iccplot.CRSM iccplot.MPRM iccplot.DRM
#' @param object Object of class \code{CRSM} for ICC of the
#' CRSM or object of class \code{MPRM} for ICC plot of the MPRM or object of class \code{DRM} for ICC plot of the DRM
#' @param items Character vector \code{"all"} to display ICC curves for all items. By entering a numeric vector, a subset of items
#'  can be chosen for which ICC plots are drawn.
#' @param \dots \dots{}
#' @author Christine Hohensinn
#' @seealso \code{\link{MPRM}} \code{\link{CRSM}} \code{\link{DRM}}
#'
#' @rdname iccplot
#' @keywords item characteristic curve, item characteristic function
#' @examples
#'
#' #estimate CRSM for the first three items
#' data(analog)
#' res_cr <- CRSM(extraversion, low=-10, high=10)
#'
#' #ICC plot
#' iccplot(res_cr)
#'
#'
#' @export iccplot
iccplot <-
  function(object,...)UseMethod("iccplot")
