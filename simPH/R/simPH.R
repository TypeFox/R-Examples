#' An R package for simulating and plotting quantities of interest from Cox
#' Proportional Hazard models.
#'
#'
#' @description An R package for simulating and plotting quantities of interest
#' (relative hazards, first differences, and hazard ratios) for linear
#' coefficients, multiplicative interactions, polynomials, penalised splines,
#' and non-proportional hazards, as well as stratified survival curves from Cox
#' Proportional Hazard models.
#'
#' The package includes the following simulation functions:
#' \itemize{
#'  \item{\code{\link{coxsimLinear}}: }{a function for simulating relative
#' hazards, first differences, hazard ratios, and hazard rates for linear,
#' non-time interacted covariates from Cox Proportional Hazard
#' (\code{\link{coxph}}) models.}
#'  \item{\code{\link{coxsimtvc}}: }{a function for simulating time-interactive
#' hazards (relative hazards, first differences, and hazard ratios) for
#' covariates from Cox Proportional Hazard models. The function will calculate
#' time-interactive hazard ratios for multiple strata estimated from a
#' stratified Cox Proportional Hazard model.}
#'  \item{\code{\link{coxsimSpline}}: }{a function for simulating quantities
#' of interest from penalised splines using multivariate normal distributions.
#' Currently does not support simulating hazard rates from stratified models.
#' Note: be extremely careful about the number of simulations you ask the
#' function to find. It is very easy to ask for more than your computer can
#' handle.}
#'  \item{\code{\link{coxsimPoly}}: }{a function for simulating quantities of
#' interest for a range of values for a polynomial nonlinear effect from Cox
#' Proportional Hazard models.}
#'  \item{\code{\link{coxsimInteract}}: }{a function for simulating quantities
#' of interest for linear multiplicative interactions, including marginal
#' effects and hazard rates.}
#' }
#' Results from these functions can be plotted using the \code{\link{simGG}}
#' method.
#'
#' The package also includes two functions that make it easier to create time
#' interactions:
#' \itemize{
#'  \item{\code{\link{SurvExpand}}: }{a function to convert a data frame of
#' non-equal interval continuous observations into equal interval continuous
#' observations.}
#'  \item{\code{\link{tvc}}: }{a function to create time interaction variables
#' that can be used in a \code{\link{coxph}} model (or any other model with time
#' interactions).}
#'  \item{setXl: }{a function for setting valid \code{Xl} values given a
#'  sequence of fitted \code{Xj} values. This makes it more intituitive to find
#'  hazard ratios and first differences for comparisons between some \eqn{Xj}
#'  fitted values and \eqn{Xl} values other than 0.}
#' }
#'
#' @references Gandrud, Christopher. 2015. simPH: An R Package for Illustrating
#' Estimates from Cox Proportional Hazard Models Including for Interactive and
#' Nonlinear Effects. Journal of Statistical Software. 65(3)1-20.
#'
#'
#' @name simPH
NULL
