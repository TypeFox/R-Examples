#' Heuristic Wright recombination for DE
#' 
#' Implements the "/wright" (Heuristic Wright) recombination for the ExpDE 
#' framework.
#'
#' @section Warning:
#' This recombination operator evaluates the candidate solutions in \code{M}, 
#' which adds an extra \code{popsize} evaluations per iteration.
#'
#' @section References:
#' F. Herrera, M. Lozano, A. M. Sanchez, "A taxonomy for the crossover
#' operator for real-coded genetic algorithms: an experimental study", 
#' International Journal of Intelligent Systems 18(3) 309-338, 2003.\cr
#' A.H. Wright, "Genetic Algorithms for Real Parameter Optimization",
#' Proc. Foundations of Genetic Algorithms, 205-218, 1991.
#'
#' @param X population matrix (original)
#' @param M population matrix (mutated) 
#' @param ... optional parameters (unused)
#' 
#' @return Matrix \code{U} containing the recombined population
#' 
#' @export

recombination_wright <- function(X, M, ...) {
  # ========== Error catching and default value definitions
  
  # Get access to variables in the calling environment
  env <- parent.frame()

  if (!identical(dim(X), dim(M))) {
    stop("recombination_wright() requires dim(X) == dim(M)")
  }
  if (!all(c("J", "probpars", "nfe") %in% names(env))){
    stop("recombination_wright() requires calling environment to contain 
         variables J, nfe and probpars")
  }
  
  # ==========
  
  # Performance values of the current population (X)
  f.X <- env$J
  
  #Evaluate population M
  f.M <- evaluate_population(probpars = env$probpars, 
                             Pop      = M)
  
  # Update NFE counter in calling environment
  env$nfe <- env$nfe + nrow(M)
  
  # Get best parent indicator matrix
  X.is.best <- matrix(rep(f.X <= f.M,
                          times = ncol(X)),
                      ncol = ncol(X),
                      byrow = FALSE)
  
  
  # Get 'best' and 'worst' parents
  C1 <- X * X.is.best + M * !X.is.best
  C2 <- M * X.is.best + X * !X.is.best
  
  # Return recombined population
  return (randM(X) * (C1 - C2) + C1)
}
