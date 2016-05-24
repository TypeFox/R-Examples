#' The posterior probability criterion function
#'
#' Thall and Simon's criterion function for determining the
#' trial decision cutoffs based on the posterior probability.
#'
#' @usage
#' postprob(y, n, alpha_e, beta_e, alpha_s, beta_s, delta)
#' @param y the number of responses among \eqn{n} patients treated by the experimental drug at a certain stage of the trial.
#' @param n the number of patients treated by the experimental drug at a certain stage of the trial.
#' @param alpha_e the hyperparameter (shape1) of the Beta prior for the experimental drug.
#' @param beta_e the hyperparameter (shape2) of the Beta prior for the experimental drug.
#' @param alpha_s the hyperparameter (shape1) of the Beta prior for the standard drug.
#' @param beta_s the hyperparameter (shape2) of the Beta prior for the standard drug.
#' @param delta the minimally acceptable increment of the response rate for the experimental drug compared with the standard drug.
#' @return
#' \item{prob}{the posterior probability: \eqn{\Pr(p_E > p_S + \delta | y)}}
#' @references
#' Thall, P. F., Simon, R. (1994).
#' Practical Bayesian guidelines for phase IIB clinical trials.
#' \emph{Biometrics} \strong{50}: 337-349.
#'
#' Yin, G. (2012).
#' \emph{Clinical Trial Design: Bayesian and Frequentist Adaptive Methods.}
#' New York: Wiley.
#' @importFrom stats integrate pbeta dbeta
#' @export
postprob <- function(y, n, alpha_e, beta_e, alpha_s, beta_s, delta) {

  f <- function(x, y, n, alpha_e, beta_e, alpha_s, beta_s, delta) {
    (1 - pbeta(x + delta, alpha_e + y, beta_e + n - y))*dbeta(x, alpha_s, beta_s)
  }

  integrate(f, 0, 1 - delta, y, n, alpha_e, beta_e, alpha_s, beta_s, delta)$value

}
