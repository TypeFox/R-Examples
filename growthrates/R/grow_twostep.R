#' Twostep Growth Model
#'
#' System of two differential equations describing bacterial growth as two-step
#' process of activation (or adaptation) and growth.
#'
#' The model is given as a system of two differential equations:
#'
#' \deqn{dy_i/dt = -kw * yi}
#' \deqn{dy_a/dt =  kw * yi + mumax * (1 - (yi + ya)/K) * ya}
#'
#' that are then numerically integrated ('simulated') according to time (t). The
#' model assumes that the population consists of active (\eqn{y_a}) and inactive
#' (\eqn{y_i}) cells so that the observed abundance is (\eqn{y = y_i + y_a}).
#' Adapting inactive cells change to the active state with a first order 'wakeup'
#' rate (\eqn{kw}).
#'
#' @param time actual time (for the ode) resp. vector of simulation time steps.
#' @param y named vector with state of the system
#'   (yi, ya: abundance of inactive and active organisms, e.g.
#'   concentration of inactive resp. active cells).
#' @param parms parameters of the two-step growth model:
#'   \itemize{
#'      \item \code{yi, ya} initial abundance of active and inactive organisms,
#'      \item \code{kw} activation (``wakeup'') constant (1/time),
#'      \item \code{mumax} maximum growth rate (1/time),
#'      \item \code{K} carrying capacity (max. abundance).
#'   }
#' @param \dots placeholder for additional parameters (for user-extended versions of this function)
#'
#' @return
#'
#' For \code{ode_twostep}: matrix containing the simulation outputs.
#' The return value of has also class \code{deSolve}.
#'
#' For \code{grow_twostep}: vector of dependent variable (\code{y}) and
#'   its log-transformed values (\code{log_y}):
#'
#' \itemize{
#' \item \code{time} time of the simulation
#' \item \code{yi} concentration of inactive cells
#' \item \code{ya} concentration of active cells
#' \item \code{y} total cell concentration
#' \item \code{log_y} natural log of total cell concentration
#' }
#'
#' @details Function \code{ode_twostep} is the system of differential equations,
#' whereas \code{grow_twostep} runs a numerical simulation over time.
#'
#' A similar two-compartment model, but without the logistic term, was discussed by Baranyi (1998).
#'
#'
#' @references
#'
#' Baranyi, J. (1998). Comparison of stochastic and deterministic concepts of bacterial lag.
#' J. heor. Biol. 192, 403--408.
#'
#' @examples
#'
#' time <- seq(0, 30, length=200)
#' parms <- c(kw = 0.1,	mumax=0.2, K=0.1)
#' y0    <-  c(yi=0.01, ya=0.0)
#' out   <- ode(y0, time, ode_twostep, parms)
#'
#' plot(out)
#'
#' o <- grow_twostep(0:100, c(yi=0.01, ya=0.0, kw = 0.1,	mumax=0.2, K=0.1))
#' plot(o)
#'
#' @family growth models
#'
#' @rdname grow_twostep
#' @export ode_twostep
#'
ode_twostep <- function (time, y, parms, ...) {
  ## the differential equations
  with(as.list(c(parms, y)), {
    dyi <- -kw * yi
    dya <-  kw * yi + mumax * (1 - (yi + ya)/K) * ya
    list(c(dyi, dya), y = unname(yi + ya), log_y = log(unname(yi + ya)))
  })
}

## @rdname grow_twostep
## @export grow_twostep.R
##
#grow_twostep.R <- function(time, parms, ...) {
#  ## assign parameters and solve differential equations
#  y0    <- parms[c("yi", "ya")]
#  parms <- parms[c("kw", "mumax", "K")]
#  out  <-  ode(y0, time, ode_twostep, parms, ...)
#}

#' @rdname grow_twostep
#' @export
#'
grow_twostep <- function(time, parms, ...) {
  ## assign parameters and solve differential equations
  #cat("compiled code running\n")
  y0    <- parms[c("yi", "ya")]
  parms <- parms[c("kw", "mumax", "K")]
  #cat(y0, "\n")
  #cat(parms, "\n")
  out <- ode(y0, time, func = "d_twostep", parms = parms,
             dllname = "growthrates",
             initfunc = "ini_twostep", nout = 2, outnames=c("y", "log_y"), ...)

}
## attach names of parameters as attributes
attr(grow_twostep, "pnames") <- c("yi","ya", "mumax", "K")
class(grow_twostep) <- c("growthmodel", "function")


