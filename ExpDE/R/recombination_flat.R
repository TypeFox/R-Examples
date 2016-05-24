#' Flat recombination for DE
#' 
#' Implements the "/flat" (Flat) recombination for the ExpDE 
#' framework
#'
#' @section References: 
#' Picek, S.; Jakobovic, D.; Golub, M., "On the recombination operator 
#' in the real-coded genetic algorithms," CEC'2013, pp.3103-3110, 2013\cr
#' F. Herrera, M. Lozano, J.L. Verdegay, "Tackling Real-Coded Genetic 
#' Algorithms: Operators and Tools for Behavioural Analysis", 
#' Artificial Intelligence Review 12 265-319, 1998.
#'
#' @param X population matrix (original)
#' @param M population matrix (mutated) 
#' @param ... optional parameters (unused)
#' 
#' @return Matrix \code{U} containing the recombined population
#' 
#' @export

recombination_flat<- function(X, M, ...) {

  # ========== Error catching and default value definitions
  if (!identical(dim(X), dim(M))) {
    stop("recombination_flat() requires dim(X) == dim(M)")
  }
  # ==========
  
  # Determine maximum and minimum values
  Cmin <- pmin(X, M)
  Cmax <- pmax(X, M)
  I    <- Cmax - Cmin
 
  # Return recombined population
  return (Cmin  + randM(X) * I)
}
