##
## control.R - Preliminary control parameter framework
##
## Author: Olaf Mersmann <olafm@statistik.tu-dortmund.de>
##

##' Basic EMOA control parameters.
##'
##' The following control parameters are recognized by \code{emoa_control}:
##' \describe{
##'   \item{logger}{\code{emoa_logger} object used to log events.}
##'   \item{n}{Number of parameters, defaults to the length of the longer
##'     of \code{upper} or \code{lower}.}
##'   \item{d}{Number of dimensions.}
##' }
##'
##' @param f Multiobjectve optimization function.
##' @param upper Upper bounds of parameter space.
##' @param lower Lower bounds of parameter space.
##' @param ... Further arguments passed to \code{f}.
##' @param control List of control parameters.
##' @param default List of default control parameters.
##'
##' @return The \code{control} list with suitably adjusted
##' arguments. Missing control parameters are taken from
##' \code{default} or, if not present there, from an internal default.
##' 
##' @export
##' @author Olaf Mersmann \email{olafm@@statistik.tu-dortmund.de}
emoa_control <- function(f, upper, lower, ..., control, default) {
  control$logger <- coalesce(control[["logger"]], default[["logger"]], emoa_console_logger())
  control$n <- as.integer(coalesce(control[["n"]],
                                   default[["n"]],
                                   max(length(lower), length(upper))))
  control$d <- as.integer(coalesce(control[["d"]],
                                   default[["d"]],
                                   length(f(rep(NA, control$n), ...))))
  control
}

##' Steady state EMOA parameters
##' 
##' \code{steady_state_emoa_control} interprets the following control
##' parameters:
##' \describe{
##'   \item{mu}{Population size.}
##'   \item{maxeval}{Maximum number of function evaluations to use.}
##' }
##'
##' @param f Multiobjectve optimization function.
##' @param upper Upper bounds of parameter space.
##' @param lower Lower bounds of parameter space.
##' @param ... Further arguments passed to \code{f}.
##' @param control List of control parameters.
##' @param default List of default control parameters.
##'
##' @return The \code{control} list with suitably adjusted
##' arguments. Missing control parameters are taken from
##' \code{default} or, if not present there, from an internal default.
##' 
##' @export
##' @author Olaf Mersmann \email{olafm@@statistik.tu-dortmund.de}
steady_state_emoa_control <- function(f, upper, lower, ..., control, default=list()) {
  control <- emoa_control(f, upper, lower, ..., control=control, default=default)
  control$mu <- as.integer(coalesce(control[["mu"]],
                                    default[["mu"]],
                                    100L))
  control$maxeval <- as.integer(coalesce(control[["maxeval"]],
                                         default[["maxeval"]],
                                         control$mu * 300L))
  control
}

##' Simulated binary crossover (SBX) control parameters
##' 
##' \code{sbx_control} interprets the following parameters used to
##' control the behaviour of the simulated binary crossover operator
##' (see \code{\link{sbx_operator}}):
##' \describe{
##'   \item{sbx.n}{Nu parameter of SBX.}
##'   \item{sbx.p}{$p$ parameter of SBX.}
##' }
##'
##' @param f Multiobjectve optimization function.
##' @param upper Upper bounds of parameter space.
##' @param lower Lower bounds of parameter space.
##' @param ... Further arguments passed to \code{f}.
##' @param control List of control parameters.
##' @param default List of default control parameters.
##'
##' @return The \code{control} list with suitably adjusted
##' arguments. Missing control parameters are taken from
##' \code{default} or, if not present there, from an internal default.
##' 
##' @export
##' @author Olaf Mersmann \email{olafm@@statistik.tu-dortmund.de}
sbx_control <- function(f, upper, lower, ..., control, default=list()) {
  if (!"crossover" %in% names(control)) {
    if (!"crossover" %in% names(default)) {
      control$sbx.n <- coalesce(control[["sbx.n"]], default[["sbx.n"]], 5)
      control$sbx.p <- coalesce(control[["sbx.p"]], default[["sbx.p"]], 1.0)
      control$crossover <- sbx_operator(control$sbx.n, control$sbx.p, lower, upper)
    } else {
      control$crossover <- default[["crossover"]]
    }
  }
  control
}

##' Polynomial muation (PM) control parameters
##'
##' Control parameters:
##' \describe{
##'   \item{pm.n}{Nu parameter of PM.}
##'   \item{pm.p}{p parameter of PM.}
##' }
##' 
##' @param f Multiobjectve optimization function.
##' @param upper Upper bounds of parameter space.
##' @param lower Lower bounds of parameter space.
##' @param ... Further arguments passed to \code{f}.
##' @param control List of control parameters.
##' @param default List of default control parameters.
##'
##' @return The \code{control} list with suitably adjusted
##' arguments. Missing control parameters are taken from
##' \code{default} or, if not present there, from an internal default.
##' 
##' @export
##' @author Olaf Mersmann \email{olafm@@statistik.tu-dortmund.de}
pm_control <- function(f, upper, lower, ..., control, default=list()) {
  if (!"mutate" %in% names(control)) {
    if (!"mutate" %in% names(default)) {
      control$pm.n <- coalesce(control[["pm.n"]], default[["pm.n"]], 10)
      control$pm.p <- coalesce(control[["pm.p"]], default[["pm.p"]], .2)
      control$mutate <- pm_operator(control$pm.n, control$pm.p, lower, upper)
    } else {
      control$mutate <- default[["mutate"]]
    }
  }
  control
}
