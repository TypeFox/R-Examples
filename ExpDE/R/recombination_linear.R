#' Linear recombination for DE
#' 
#' Implements the "/linear" recombination for the ExpDE framework
#'
#' @section Warning:
#' This recombination operator evaluates \code{3*popsize} candidate solutions 
#' per iteration of the algorithm. The value of the \code{nfe} counter and the 
#' vector of performance values \code{G} are updated in the calling environment.
#' 
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

recombination_linear <- function(X, M, ...) {
  # ========== Error catching and default value definitions
  
  # Get access to variables in the calling environment
  env <- parent.frame()
  
  if (!identical(dim(X), dim(M))) {
    stop("recombination_linear() requires dim(X) == dim(M)")
  }
  if (!all(c("J", "probpars", "nfe") %in% names(env))){
    stop("recombination_linear() requires calling environment to contain 
         variables J, nfe and probpars")
  }
  
  # ==========
  # Generate trial offspring 
  H1 <- (0.5 * X) + (0.5 * M)
  H2 <- (1.5 * X) - (0.5 * M)
  H3 <- -(0.5 * X) + (1.5 * M)
  
  # Evaluate trial offspring
  f1 <- evaluate_population(probpars = env$probpars, 
                            Pop      = H1)
  
  f2 <- evaluate_population(probpars = env$probpars, 
                            Pop      = H2)
  
  f3 <- evaluate_population(probpars = env$probpars, 
                            Pop      = H3)
  
  env$nfe <- env$nfe + 3 * nrow(X)
  
  # Perform recombination
  fbest <- pmin(f1, f2, f3)
  
  # Update performance vector in calling environment
  env$G[f1 == fbest] <- f1[f1 == fbest]
  env$G[f2 == fbest] <- f2[f2 == fbest]
  env$G[f3 == fbest] <- f3[f3 == fbest]
  
  Pop.trialx <- X
  Pop.trialx[f1==fbest, ] <- H1[f1==fbest, ]
  Pop.trialx[f2==fbest, ] <- H2[f2==fbest, ]
  Pop.trialx[f3==fbest, ] <- H3[f3==fbest, ]
  
  # Return recombined population
  return (Pop.trialx)
  }
