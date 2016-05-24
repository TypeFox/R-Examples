#' The Baranyi and Roberts Growth Model
#'
#' The growth model of Baranyi and Roberts (1995) written as analytical solution
#' of the system of differential equations.
#'
#' The version of the equation used in this package has the following form:
#'
#' \deqn{A = time + 1/mumax * log(exp(-mumax * time) + exp(-h0) - exp(-mumax * time - h0))}
#' \deqn{log(y) = log(y0) + mumax * A - log(1 + (exp(mumax * A) - 1) / exp(log(K) - log(y0)))}
#'
#' @param time vector of time steps (independent variable).
#' @param parms named parameter vector of the Baranyi growth model with:
#' \itemize{
#'   \item \code{y0} initial value of abundance,
#'   \item \code{mumax} maximum growth rate (1/time),
#'   \item \code{K} carrying capacity (max. abundance),
#'   \item \code{h0} parameter specifying the initial physiological state of
#'     organisms (e.g. cells) and in consequence the lag phase
#'     (h0 = max growth rate * lag phase).
#' }
#'
#' @return vector of dependent variable (\code{y}) and its log-transformed
#'   values (\code{log_y}).
#'
#'
#'
#' @references
#'
#' Baranyi, J. and Roberts, T. A. (1994).
#' A dynamic approach to predicting bacterial growth in food.
#' International Journal of Food Microbiology, 23, 277-294.
#'
#' Baranyi, J. and Roberts, T.A. (1995). Mathematics of predictive microbiology.
#' International Journal of Food Microbiology, 26, 199-218.
#'
#'
#' @examples
#'
#' time <- seq(0, 30, length=200)
#' y    <- grow_baranyi(time, c(y0=0.01, mumax=.5, K=0.1, h0=5))[,"y"]
#' plot(time, y, type="l")
#' plot(time, y, type="l", log="y")
#'
#' @family growth models
#'
#' @export
#'
grow_baranyi <- function(time, parms) {
  with(as.list(parms), {
    ## todo: q0 in original paper, h0 in Huang
    A <- time + 1/mumax * log(exp(-mumax * time) + exp(-h0) - exp(-mumax * time - h0))
    #log_y <- y0 + mumax * A - log(1 + (exp(mumax * A) - 1)/(exp(K - y0)))
    log_y <- log(y0) + mumax * A - log(1 + (exp(mumax * A) - 1) / exp(log(K) - log(y0)))

    return(as.matrix(data.frame(time = time, y = exp(log_y), log_y = log_y)))
  })
}
## attach names of parameters as attributes
attr(grow_baranyi, "pnames") <- c("y0", "mumax", "K", "h0")
class(grow_baranyi) <- c("growthmodel", "function")

