#' @title Log-likelihood of LICORS model
#'
#' @description 
#' Computes the \emph{average} log-likelihood \eqn{\frac{1}{N} \ell(\mathbf{W}; \mathcal{D})} as a function of the 
#' weight matrix \eqn{\mathbf{W}} and the predictive state distributions 
#' \eqn{P(X = x \mid S = s_j) \approx P(X = x \mid \mathbf{W}_j)} for all 
#' \eqn{j = 1, \ldots, K}. See References.
#'
#' @param weight.matrix \eqn{N \times K} weight matrix
#' @param pdfs.FLC an \eqn{N \times K} matrix containing the estimates of all \eqn{K} 
#' FLC densities evaluated at all \eqn{N} sample FLCs.
#' @param lambda regularization parameter. Default: \code{lambda=0} 
#' (\code{penalty} and \code{q} will be ignored in this case).
#' @param penalty type of penalty: \code{c("entropy", "1-Lq", "lognorm")}. 
#' Default: \code{"entropy"}
#' @param base logarithm base for the \code{"entropy"} penalty.
#' Default: \code{base = 2}.  Any other real number is allowed; 
#' if \code{base = "num.states"} then it will internally assign it 
#' \code{base = ncol(weight.matrix)}.
#' @param q exponent for \eqn{L_q} norm.
#' @keywords manip nonparametric
#' @export
#' 

compute_LICORS_loglik <- function(weight.matrix, pdfs.FLC, 
                                  lambda = 0, penalty = "entropy", q = 2, 
                                  base = exp(1)) {
  qq <- q
  loglik <- mean(log(rowSums(pdfs.FLC * weight.matrix)))
  if (lambda > 0) {
    #NN <- nrow(pdfs.FLC)
    #kk <- ncol(pdfs.FLC)
    loglik <- loglik - lambda * compute_mixture_penalty(weight.matrix, 
                                                        type = penalty, 
                                                        q = qq, base = base) #- 1/2 * kk * log(NN)/NN
  }
  return(loglik)
} 