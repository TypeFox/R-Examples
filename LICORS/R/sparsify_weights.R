#' @title Sparsify weights
#'
#' @description 
#' This function makes weights of a mixture model more sparse using gradient 
#' based penalty methods.
#' 
#' @param weight.matrix.proposed \eqn{N \times K} weight matrix
#' @param weight.matrix.current \eqn{N \times K} weight matrix
#' @param penalty type of penalty: \code{c("entropy", "1-Lq", "lognorm")}. 
#' Default: \code{"entropy"}
#' @param lambda penalization parameter: larger \code{lambda} gives sparser 
#' mixture weights
#' @keywords manip array
#' @export
#' @seealso \code{\link{compute_mixture_penalty}}, \code{\link{mixed_LICORS}}
#' @examples
#' WW = matrix(c(rexp(10, 1/10), runif(10)), ncol=5, byrow = FALSE)
#' WW = normalize(WW)
#' WW_sparse = sparsify_weights(WW, lambda = 0.1)
#' WW_more_sparse = sparsify_weights(WW, lambda = 0.5)
#' compute_mixture_penalty(WW)
#' compute_mixture_penalty(WW_sparse)
#' compute_mixture_penalty(WW_more_sparse)

sparsify_weights <- function(weight.matrix.proposed, 
                             weight.matrix.current = NULL, 
                             penalty = "entropy", 
                             lambda = 0) {
  
  if (lambda == 0) {
    return(weight.matrix.proposed)
  } else {
    kk <- ncol(weight.matrix.proposed)
    if (is.null(weight.matrix.current)) {
      penalty.gradient <- - (1 + log2(weight.matrix.proposed + 10^(-3))) / log2(kk)
    } else {
      penalty.gradient <- - (1 + log2(weight.matrix.current + 10^(-3))) / log2(kk)
    }
    # lambda = lambda + runif(1, -lambda/10, lambda/10)
    
    # without hessian weight.matrix.new = weight.matrix.proposed + lambda *
    # penalty.gradient multiply by Hessian
    weight.matrix.new <- weight.matrix.proposed * (1 - log2(kk) * lambda * penalty.gradient)
    weight.matrix.new[weight.matrix.new < 0] <- 0
    all.zeros <- apply(weight.matrix.new, 1, function(u) all(u == 0))
    
    if (sum(all.zeros) > 1) {
      weight.matrix.new[all.zeros, apply(weight.matrix.proposed[all.zeros, ], 1, which.max)] <- 1
    } else {
      weight.matrix.new[all.zeros, which.max(weight.matrix.proposed[all.zeros, ])] <- 1
    }
    return(normalize(weight.matrix.new))
  }
} 