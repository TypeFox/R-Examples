#' Exponential Growth Model
#'
#' Unlimited exponential growth model.
#'
#' The equation used is:
#' \deqn{y = y0 * exp(mumax * time)}
#'
#' @param time vector of time steps (independent variable).
#' @param parms named parameter vector of the exponential growth model with:
#' \itemize{
#'   \item \code{y0} initial abundance (e.g. concentration of bacterial cells).
#'   \item \code{mumax} maximum growth rate (1/time).
#' }
#'
#' @return vector of dependent variable (\code{y}) and its log-transformed
#'   values (\code{log_y}).
#'
#' @examples
#'
#' time <- seq(0, 30, length=200)
#' y <- grow_exponential(time, c(y0=1, mumax=0.5))[,"y"]
#' plot(time, y, type="l")
#'
#' @family growth models
#'
#' @rdname grow_exponential
#' @export
#'
grow_exponential <- function(time, parms) {
  ## lm object coefficients have no names
  y0 <- parms[1]
  mumax <- parms[2]
  y  <- y0 * exp(mumax * time)
  return(as.matrix(data.frame(time=time, y=y, log_y=log(y))))
}
## attach names of parameters as attributes
attr(grow_exponential, "pnames") <- c("y0", "mumax")
class(grow_exponential) <- c("growthmodel", "function")
