#' @title Compute Negative Entropy Criterion (NEC)
#'
#' @description 
#' Computes the negative entropy criterion (NEC) to assess the number of clusters
#' in a mixture model. See References for details.
#'
#' @param weight.matrix \eqn{N \times K} weight matrix
#' @param loglik.1 baseline log-likelihood for \eqn{K=1} cluster model
#' @param loglik.k log-likelihood for \eqn{K} cluster model
#' @references
#' Christophe Biernacki, Gilles Celeux, and G\'erand Govaert(1999). 
#' ``An improvement of the NEC criterion for assessing the number of clusters 
#' in a mixture model''. Non-Linear Anal. 20, 3, 267-272.
#' @keywords manip array
#' @export
#' @seealso \code{\link{compute_mixture_penalty}}
#' @examples
#' WW = matrix(c(rexp(10, 1/10), runif(10)), ncol = 5, byrow=FALSE)
#' WW = normalize(WW)
#' compute_NEC(WW, -2, -1)

compute_NEC <- function(weight.matrix, loglik.1 = NULL, loglik.k = NULL){
  if (is.null(loglik.1) || is.null(loglik.k)){
    stop("You must provide the log-likelihoods for 1-cluster and (!) 
         K-cluster solution.")
  }
  diff.in.loglik <- loglik.k - loglik.1
  nec <- compute_mixture_penalty(weight.matrix, "entropy") / diff.in.loglik
  return(nec)
}
