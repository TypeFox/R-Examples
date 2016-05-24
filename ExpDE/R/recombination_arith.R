#' Arithmetic recombination for DE
#' 
#' Implements the "/arith" (arithmetic) recombination for the ExpDE framework
#' 
#' @section References:
#' F. Herrera, M. Lozano, A. M. Sanchez, "A taxonomy for the crossover
#' operator for real-coded genetic algorithms: an experimental study", 
#' International Journal of Intelligent Systems 18(3) 309-338, 2003.
#'
#' @param X population matrix (original)
#' @param M population matrix (mutated) 
#' @param ... optional parameters (unused)
#' 
#' @return Matrix \code{U} containing the recombined population
#' 
#' @export

recombination_arith <- function(X, M, ...) {
  
  # ========== Error catching and default value definitions
  if (!identical(dim(X), dim(M))) {
    stop("recombination_arith() requires dim(X) == dim(M)")
  }
  # ==========
  
  lambda <- matrix(rep(stats::runif(nrow(X)),
                       ncol(X)),
                   ncol = ncol(X),
                   byrow = FALSE)
  
  # Return recombined population
  return(X*lambda + M*(1 - lambda))
}
