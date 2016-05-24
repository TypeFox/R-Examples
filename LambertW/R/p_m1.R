#' @title Non-principal branch probability
#' 
#' @description Computes the probability that (at least) one (out of n)
#'     observation(s) of the latent variable \eqn{U} lies in the non-principal
#'     branch region. The '\code{m1}' in \code{p_m1} stands for 'minus 1', i.e,
#'     the non-principal branch.
#' 
#' See Goerg (2011) and Details for mathematical derivations.
#' 
#' @details
#' The probability that one observation of the latent RV U lies in the
#' non-principal region equals at most \deqn{ p_{-1}(\gamma, n=1)
#' = P\left(U < -\frac{1}{|\gamma|}\right), } where \eqn{U} is the zero-mean,
#' unit variance version of the input \eqn{X \sim F_X(x \mid \boldsymbol
#' \beta)} -- see References.
#' 
#' For \eqn{N} independent RVs \eqn{U_1, \ldots, U_N}, the probability that at
#' least one data point came from the non-principal region equals
#' \deqn{
#' p_{-1}(\gamma, n=N) = P\left(U_i < -\frac{1}{|\gamma|} \; for \; at \;
#' least \; one \; i \right) } 
#' 
#' This equals (assuming independence) 
#' \deqn{ P\left(U_i < -\frac{1}{|\gamma|} \; for \; at
#' \; least \; one \; i \right) = 1 - P\left(U_i \geq -\frac{1}{|\gamma|},
#' \forall i \right) = 1 - \prod_{i=1}^{N} P\left(U_i \geq -\frac{1}{|\gamma|}
#' \right) } \deqn{ = 1 - \prod_{i=1}^{N} \left(1 - p_{-1}(\gamma, n=1) \right)
#' = 1 - (1-p_{-1}(\gamma, n=1))^N. }
#' 
#' For improved numerical stability the cdf of a geometric RV
#' (\code{\link[stats]{pgeom}}) is used to evaluate the last
#' expression. Nevertheless, numerical problems can occur for \eqn{|\gamma| <
#' 0.03} (returns \code{0} due to rounding errors).
#' 
#' Note that \eqn{1 - (1-p_{-1}(\gamma, n=1))^N} reduces to \eqn{p_{-1}(\gamma)}
#' for \eqn{N=1}.
#' 
#' @inheritParams common-arguments
#' @param gamma scalar; skewness parameter.
#' @param n number of RVs/observations.
#' @return 
#' non-negative float; the probability \eqn{p_{-1}} for \code{n} observations.
#' @keywords univar
#' @export
#' @examples
#' 
#' beta.01 <- c(mu = 0, sigma = 1)
#' # for n=1 observation
#' p_m1(0, beta = beta.01, distname = "normal") # identical to 0
#' # in theory != 0; but machine precision too low
#' p_m1(0.01, beta = beta.01, distname = "normal") 
#' p_m1(0.05, beta = beta.01, distname = "normal") # extremely small
#' p_m1(0.1, beta = beta.01, distname = "normal") # != 0, but very small
#' # 1 out of 4 samples is a non-principal input;
#' p_m1(1.5, beta = beta.01, distname = "normal") 
#' # however, gamma=1.5 is not common in practice
#' 
#' # for n=100 observations
#' p_m1(0, n=100, beta = beta.01, distname = "normal") # == 0
#' p_m1(0.1, n=100, beta = beta.01, distname = "normal") # still small
#' p_m1(0.3, n=100, beta = beta.01, distname = "normal") # a bit more likely
#' p_m1(1.5, n=100, beta = beta.01, distname = "normal") 
#' # Here we can be almost 100% sure (rounding errors) that at least one
#' # y_i was caused by an input in the non-principal branch.
#'

p_m1 <- function(gamma, beta, distname, n = 1, use.mean.variance = TRUE) {
  
  stopifnot(n > 0,
            n == round(n),
            is.numeric(gamma),
            length(gamma) == 1)
  check_distname(distname)
  check_beta(beta, distname = distname)
  
  if (gamma < 0) {
    gamma <- -gamma
  }
  out <- pU(-1/gamma, beta = beta, distname = distname,
            use.mean.variance = use.mean.variance)
  names(out) <- NULL
  
  if (n > 1 && !is.na(out)) {
    out <- suppressWarnings(pgeom(n - 1, prob = out))
    if (is.na(out)) {
      out <- 0
    }
  }
  return(out)
} 
