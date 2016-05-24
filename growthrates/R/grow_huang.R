#' Growth Model According to Huang
#'
#' Huangs growth model written as analytical solution of the differential equations.
#'
#' The version of the equation used in this package has the following form:
#' \deqn{B = time + 1/alpha * log((1+exp(-alpha * (time - lambda)))/(1 + exp(alpha * lambda)))}
#' \deqn{log(y) = log(y0) + log(K) - log(y0 + (K - y0) * exp(-mumax * B))}
#'
#' In contrast to the original publication, all parameters related to population
#' abundance (y, y0, K) are given as untransformed values.
#' They are not log-transformed.\cr
#' In general, using log-transformed parameters would indeed be a good idea to
#' avoid the need of constained optimization, but tests showed that
#' box-constrained optimization worked resonably well.
#' Therefore, handling of optionally log-transformed parameters was removed
#' from the package to avoid confusion. If you want to discuss this, please
#' let me know.
#'
#' @param time vector of time steps (independent variable).
#' @param parms named parameter vector of Huang's growth model with:
#' \itemize{
#'   \item \code{y0} initial value of abundance,
#'   \item \code{mumax} maximum growth rate (1/time),
#'   \item \code{K} carrying capacity (max. total concentration of cells),
#'   \item \code{alpha} shape parameter determining the curvature,
#'   \item \code{lambda} parameter determining the lag time.
#'
#' }
#'
#' @return vector of dependent variable (\code{y}) and its log-transformed
#'   values (\code{log_y}).
#'
#'
#' @references
#'
#' Huang, Lihan (2008) Growth kinetics of Listeria monocytogenes in broth and
#' beef frankfurters - determination of lag phase duration and exponential
#' growth rate under isothermal conditions. Journal of Food Science 73(5),
#' E235 -- E242. doi:10.1111/j.1750-3841.2008.00785.x
#'
#' Huang, Lihan (2011) A new mechanistic growth model for simultaneous
#' determination of lag phase duration and exponential growth rate and a new
#' Belehdradek-type model for evaluating the effect of temperature on growth rate.
#' Food Microbiology 28, 770 -- 776. doi:10.1016/j.fm.2010.05.019
#'
#' Huang, Lihan (2013) Introduction to USDA Integrated Pathogen Modeling
#' Program (IPMP). Residue Chemistry and Predictive Microbiology Research
#' Unit. USDA Agricultural Research Service.
#'
#'
#'
#' @examples
#'
#' time <- seq(0, 30, length=200)
#' y    <- grow_huang(time, c(y0=0.01, mumax=.1, K=0.1, alpha=1.5, lambda=3))[,"y"]
#' plot(time, y, type="l")
#' plot(time, y, type="l", log="y")
#'
#' @family growth models
#'
#' @rdname grow_huang
#' @export
#'
grow_huang <- function(time, parms) {
  with(as.list(parms), {
    B <- time + 1/alpha * log((1+exp(-alpha * (time - lambda)))/(1 + exp(alpha * lambda)))
    #log_y <- y0 + K - log(exp(y0) + (exp(K) - exp(y0)) * exp(-mumax * B))
    log_y <- log(y0) + log(K) - log(y0 + (K - y0) * exp(-mumax * B))
    return(as.matrix(data.frame(time = time, y = exp(log_y), log_y = log_y)))
  })
}
## attach names of parameters as attributes
attr(grow_huang, "pnames") <- c("y0", "mumax", "K", "alpha", "lambda")
class(grow_huang) <- c("growthmodel", "function")
