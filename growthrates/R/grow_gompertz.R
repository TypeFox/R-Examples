#' Growth Model According to Gompertz
#'
#' Gompertz growth model written as analytical solution of the differential
#'   equation system.
#'
#' The equation used here is:
#' \deqn{y = K * exp(log(y0 / K) * exp(-mumax * time))}
#'
#' @param time vector of time steps (independent variable).
#' @param parms named parameter vector of the Gompertz growth model with:
#' \itemize{
#'   \item \code{y0} initial value of abundance,
#'   \item \code{mumax} maximum growth rate (1/time),
#'   \item \code{K} maximum abundance (carrying capacity).
#'
#' }
#'
#'
#' @references Tsoularis, A. (2001) Analysis of Logistic Growth Models.
#'   Res. Lett. Inf. Math. Sci, (2001) 2, 23-46.
#'
#' @return vector of dependent variable (\code{y})
#'
#' @examples
#'
#' time <- seq(0, 30, length=200)
#' y    <- grow_gompertz(time, c(y0=1, mumax=.2, K=10))[,"y"]
#' plot(time, y, type="l", ylim=c(0, 20))
#'
#'
#' @family growth models
#'
#' @rdname grow_gompertz
#' @export grow_gompertz
#'
grow_gompertz <- function(time, parms) {
  with(as.list(parms), {
    y <- K * exp(log(y0 / K) * exp(-mumax * time))
    return(as.matrix(data.frame(time = time, y = y, log_y = log(y))))
  })
}
## attach names of parameters as attributes
attr(grow_gompertz, "pnames") <- c("y0", "mumax", "K")
class(grow_gompertz) <- c("growthmodel", "function")
