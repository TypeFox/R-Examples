#' Growth Model According to Richards
#'
#' Richards growth model written as analytical solution of the differential equation.
#'
#' The equation used is:
#'
#' \deqn{y = K*(1-exp(-beta * mumax * time)*(1-(y0/K)^-beta))^(-1/beta)}
#'
#' The naming of parameters used here follows the convention of Tsoularis (2001),
#' but uses \code{mumax} for growtrate and \code{y} for abundance to make them
#' consistent to other growth functions.
#'
#' @param time vector of time steps (independent variable).
#' @param parms named parameter vector of the Richards growth model with:
#' \itemize{
#'   \item \code{y0} initial value of abundance,
#'   \item \code{mumax} maximum growth rate (note different interpretation compared
#'     to exponential growth),
#'   \item \code{K} carrying capacity (max. total concentration of cells),
#'   \item \code{beta} shape parameter determining the curvature.
#'
#' }
#'
#' @return vector of dependent variable (\code{y}) and its log-transformed
#'   values (\code{log_y}).
#'
#' @references
#'
#' Richards, F. J. (1959) A Flexible Growth Function for Empirical Use.
#' Journal of Experimental Botany 10 (2): 290--300.
#'
#' Tsoularis, A. (2001) Analysis of Logistic Growth Models.
#' Res. Lett. Inf. Math. Sci, (2001) 2, 23--46.
#'
#' @examples
#'
#' time <- seq(0, 30, length=200)
#' y    <- grow_richards(time, c(y0=1, mumax=.5, K=10, beta=2))[,"y"]
#' plot(time, y, type="l")
#' y    <- grow_richards(time, c(y0=1, mumax=.5, K=10, beta=100))[,"y"]
#' lines(time, y, col="red")
#' y    <- grow_richards(time, c(y0=1, mumax=.5, K=10, beta=.2))[,"y"]
#' lines(time, y, col="blue")
#'
#' @family growth models
#'
#' @rdname grow_richards
#' @export grow_richards
#'
grow_richards <- function(time, parms) {
  with(as.list(parms), {
    y <- K*(1-exp(-beta * mumax * time)*(1-(y0/K)^-beta))^(-1/beta)

    return(as.matrix(data.frame(time = time, y = y, log_y = log(y))))
  })
}
## attach names of parameters as attributes
attr(grow_richards, "pnames") <- c("y0", "mumax", "K", "beta")
class(grow_richards) <- c("growthmodel", "function")
