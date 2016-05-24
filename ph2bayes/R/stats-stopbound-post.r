#' The stopping boundaries based on Thall and Simon's criterion
#'
#' The stopping boundaries based on Thall and Simon's criterion.
#'
#' @usage
#' stopbound_post(theta, type, nmax, alpha_e, beta_e, alpha_s, beta_s, delta)
#' @param theta the cutoff probability: typically, \eqn{\theta = [0.95, 0.99]} for superiority, \eqn{\theta = [0.01, 0.05]} for futility.
#' @param type type of boundaries: "superiority" or "futility".
#' @param nmax the maximum number of patients treated by the experimental drug.
#' @param alpha_e the hyperparameter (shape1) of the Beta prior for the experimental drug.
#' @param beta_e the hyperparameter (shape2) of the Beta prior for the experimental drug.
#' @param alpha_s the hyperparameter (shape1) of the Beta prior for the standard drug.
#' @param beta_s the hyperparameter (shape2) of the Beta prior for the standard drug.
#' @param delta the minimally acceptable increment of the response rate for the experimental drug compared with the standard drug.
#' Note: if type = "superiority", then delta is set to 0.
#' @return
#' \item{boundset}{the boundaries set; \eqn{U_n} or \eqn{L_n}}
#' @references
#' Thall, P. F., Simon, R. (1994).
#' Practical Bayesian guidelines for phase IIB clinical trials.
#' \emph{Biometrics} \strong{50}: 337-349.
#'
#' Yin, G. (2012).
#' \emph{Clinical Trial Design: Bayesian and Frequentist Adaptive Methods.}
#' New York: Wiley.
#' @examples
#' stopbound_post(0.05, "futility", 40, 0.6, 1.4, 15, 35, 0)
#' stopbound_post(0.05, "futility", 30, 0.4, 1.6, 10, 40, 0)
#' stopbound_post(0.95, "superiority", 40, 0.6, 1.4, 15, 35, 0)
#' @export
stopbound_post <- function(theta, type = c("superiority", "futility"), nmax, alpha_e, beta_e, alpha_s, beta_s, delta) {

  type <- match.arg(type)

  if (type == "superiority") {
    delta <- 0
  }

  bound <- rep(NA, nmax)
  for (n in 1:nmax) {
    if (type == "superiority") {
      for (y in seq(0, n)) {
        if (postprob(y, n, alpha_e, beta_e, alpha_s, beta_s, delta) >= theta) {
          bound[n] <- y
          break
        }
      }
    } else {
      for (y in seq(n, 0)) {
        if (postprob(y, n, alpha_e, beta_e, alpha_s, beta_s, delta) <= theta) {
          bound[n] <- y
          break
        }
      }
    }
  }

  boundset <- data.frame(n = 1:nmax, bound = bound)

  return(boundset[!duplicated(boundset[, 2]), ])

}

