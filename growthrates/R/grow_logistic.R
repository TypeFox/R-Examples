#' Logistic Growth Model
#'
#' Classical logistic growth model written as analytical solution of the differential equation.
#'
#' The equation used is:
#' \deqn{y = (K * y0) / (y0 + (K - y0) * exp(-mumax * time))}
#'
#' @param time vector of time steps (independent variable)
#' @param parms named parameter vector of the logistic growth model with:
#' \itemize{
#'   \item \code{y0} initial value of population measure
#'   \item \code{mumax} intrinsic growth rate (1/time)
#'   \item \code{K} carrying capacity (max. total concentration of cells)
#' }
#'
#' @return vector of dependent variable (\code{y}) and its log-transformed
#'   values (\code{log_y}).
#'
#' @examples
#'
#' time <- seq(0, 30, length=200)
#' y    <- grow_logistic(time, c(y0=1, mumax=0.5, K=10))[,"y"]
#' plot(time, y, type="l")
#'
#' @family growth models
#'
#' @rdname grow_logistic
#' @export grow_logistic
#'
grow_logistic <- function(time, parms) {
  with(as.list(parms), {
    y <- (K * y0) / (y0 + (K - y0) * exp(-mumax * time))
    return(as.matrix(data.frame(time=time, y=y, log_y=log(y))))
  })
}
## attach names of parameters as attributes
attr(grow_logistic, "pnames") <- c("y0", "mumax", "K")
class(grow_logistic) <- c("growthmodel", "function")
