#' Generalized Logistic Growth Model
#'
#' Generalized logistic growth model solved as differential equation.
#'
#' The model is given as its first derivative:
#'
#' \deqn{dy/dt = mumax * y^alpha * (1-(y/K)^beta)^gamma}
#'
#' that is then numerically integrated ('simulated') according to time (t).
#'
#' @param time vector of simulation time steps
#' @param y named vector with initial value of the system (e.g. cell concentration)
#' @param parms parameters of the generalized logistic growth model
#'   \itemize{
#'      \item \code{mumax} maximum growth rate (1/time)
#'      \item \code{K} carrying capacity (max. abundance)
#'      \item \code{alpha, beta, gamma} parameters determining the shape of growth.
#'        Setting all values to one returns the ordinary logistic function.
#'   }
#' @param \dots additional parameters passed to the \code{ode}-function.
#'
#' @return
#'
#' For \code{ode_genlogistic}: matrix containing the simulation outputs.
#' The return value of has also class \code{deSolve}.
#'
#' For \code{grow_genlogistic}: vector of dependent variable (\code{y}) and
#'   its log-transformed values (\code{log_y}).
#'
#' \itemize{
#' \item \code{time} time of the simulation
#' \item \code{y} abundance of organisms
#' \item \code{log_y} natural log of abundance
#' }
#'
#' @details The generalized logistic according to Tsoularis (2001) is a flexible
#'   model that covers exponential and logistic growth, Richards, Gompertz, von
#'   Bertalanffy, and some more as special cases.
#'
#'   The differential equation is solved numerically, where function
#'   \code{ode_genlogistic} is the differential equation, and
#'   \code{grow_genlogistic} runs a numerical simulation over time.
#'
#'   The default version \code{grow_genlogistic} is run directly as compiled code,
#'   whereas the R versions \code{ode_logistic} is
#'   provided for testing by the user.
#'
#' @references
#' Tsoularis, A. (2001) Analysis of Logistic Growth Models.
#' Res. Lett. Inf. Math. Sci, (2001) 2, 23-46.
#'
#' @examples
#'
#' time <- seq(0, 30, length=200)
#' parms <- c(mumax=0.5, K=10, alpha=1, beta=1, gamma=1)
#' y0    <-  c(y=.1)
#' out   <- ode(y0, time, ode_genlogistic, parms)
#' plot(out)
#'
#' out2 <- ode(y0, time, ode_genlogistic, parms = c(mumax=0.2, K=10, alpha=2, beta=1, gamma=1))
#' out3 <- ode(y0, time, ode_genlogistic, parms = c(mumax=0.2, K=10, alpha=1, beta=2, gamma=1))
#' out4 <- ode(y0, time, ode_genlogistic, parms = c(mumax=0.2, K=10, alpha=1, beta=1, gamma=2))
#' out5 <- ode(y0, time, ode_genlogistic, parms = c(mumax=0.2, K=10, alpha=.5, beta=1, gamma=1))
#' out6 <- ode(y0, time, ode_genlogistic, parms = c(mumax=0.2, K=10, alpha=1, beta=.5, gamma=1))
#' out7 <- ode(y0, time, ode_genlogistic, parms = c(mumax=0.3, K=10, alpha=1, beta=1, gamma=.5))
#' plot(out, out2, out3, out4, out5, out6, out7)
#'
#' ## growth with lag (cf. log_y)
#' plot(ode(y0, time, ode_genlogistic, parms = c(mumax=1, K=10, alpha=2, beta=.8, gamma=5)))
#'
#'
#' @family growth models
#'
#' @rdname grow_genlogistic
#' @export ode_genlogistic
#'
ode_genlogistic <- function (time, y, parms, ...) {
  ## the differential equations
  with(as.list(c(parms)), {
    dy <-  mumax * y^alpha * (1-(y/K)^beta)^gamma
    list(dy, log_y = log(unname(y)))
  })
}

# @rdname grow_genlogistic
# @export grow_genlogistic.R
#
# grow_genlogistic.R <- function(time, parms, ...) {
#   ## assign parameters and solve differential equations
#   y0    <- c(y = unname(parms[c("y0")]))
#   parms <- parms[c("mumax", "K", "alpha", "beta", "gamma")]
#   out  <-  as.matrix(ode(y0, time, ode_genlogistic, parms, ...))
# }

#' @rdname grow_genlogistic
#' @export grow_genlogistic
#'
grow_genlogistic <- function(time, parms, ...) {
  ## assign parameters and solve differential equations
  #cat("compiled code running\n")
  y0    <- c(y = unname(parms[c("y0")]))
  parms <- parms[c("mumax", "K", "alpha", "beta", "gamma")]
  #cat(y0, "\n")
  #cat(parms, "\n")
  out <- ode(y0, time, func = "d_genlogistic", parms = parms,
             dllname = "growthrates",
             initfunc = "ini_genlogistic", nout = 0, ...)
  cbind(out, log_y = log(out[,"y"]))
}
## attach names of parameters as attributes
attr(grow_genlogistic, "pnames") <- c("y0", "mumax", "K", "alpha", "beta", "gamma")
class(grow_genlogistic) <- c("growthmodel", "function")


