#' Cost Function for Nonlinear Fits
#'
#' Defines a cost function from the residual sum of squares between model
#' and observational data.
#'
#'
#' @param FUN function of growth model to be fitted.
#' @param p vector of fitted parameters of the growth model.
#' @param fixed.p vector of fixed  parameters of the growth model.
#' @param \dots additional parameters passed to the optimizer.
#'
#' @return an object of class \code{modCost}, see \code{\link{modCost}} in
#'   package \pkg{FME}
#'
#' @export cost
#' @keywords internal
#'
#' @details
#'
#' Function 'cost' is implemented as follows, see package FME for details:
#' \preformatted{
#' cost <- function(p, obs, FUN, fixed.p = NULL, ...) {
#'   out <- FUN(obs$time, c(p, fixed.p))
#'   modCost(out, obs, weight = "none", ...)
#' }
#' }

cost <- function(p, obs, FUN, fixed.p = NULL, ...) {
  out <- FUN(obs$time, c(p, fixed.p))
  modCost(out, obs, weight = "none", ...)
}
