#' The stopping boundaries based on Lee and Liu's criterion
#'
#' The stopping boundaries based on Lee and Liu's criterion.
#'
#' @usage
#' stopbound_pred(theta, type, nmax, alpha_e, beta_e, p_s, theta_t)
#' @param theta the cutoff probability: typically, \eqn{\theta = [0.95, 0.99]} for superiority, \eqn{\theta = [0.01, 0.05]} for futility.
#' @param type type of boundaries: "superiority" or "futility".
#' @param nmax the maximum number of patients treated by the experimental drug.
#' @param alpha_e the hyperparameter (shape1) of the Beta prior for the experimental drug.
#' @param beta_e the hyperparameter (shape2) of the Beta prior for the experimental drug.
#' @param p_s the the response rate for the standard drug.
#' @param theta_t the prespecified target probability; tipically, \eqn{\theta_T = [0.85, 0.95]}.
#' @return
#' \item{boundset}{the boundaries set: \eqn{U_n} or \eqn{L_n}}
#' @references
#' Lee, J. J., Liu, D. D. (2008).
#' A predictive probability design for phase II cancer clinical trials.
#' \emph{Clinical Trials} \strong{5}: 93-106.
#'
#' Yin, G. (2012).
#' \emph{Clinical Trial Design: Bayesian and Frequentist Adaptive Methods.}
#' New York: Wiley.
#' @examples
#' stopbound_pred(0.05, "futility", 40, 0.6, 1.4, 0.3, 0.85)
#' stopbound_pred(0.05, "futility", 30, 0.4, 1.6, 0.2, 0.85)
#' stopbound_pred(0.95, "superiority", 40, 0.6, 1.4, 0.3, 0.85)
#' @export
stopbound_pred <- function(theta, type = c("superiority", "futility"), nmax, alpha_e, beta_e, p_s, theta_t) {

  type <- match.arg(type)

  bound <- rep(NA, nmax)
  for (n in 1:nmax) {
    if (type == "superiority") {
      for (y in seq(0, n)) {
        if (predprob(y, n, nmax, alpha_e, beta_e, p_s, theta_t) >= theta) {
          bound[n] <- y
          break
        }
      }
    } else {
      for (y in seq(n, 0)) {
        if (predprob(y, n, nmax, alpha_e, beta_e, p_s, theta_t) <= theta) {
          bound[n] <- y
          break
        }
      }
    }
  }

  boundset <- data.frame(n = 1:nmax, bound = bound)

  return(boundset[!duplicated(boundset[, 2]), ])

}
