#' Blend Alpha Beta recombination for DE
#' 
#' Implements the "/blxAlphaBeta" (Blend Alpha Beta) recombination for the ExpDE 
#' framework
#' 
#' @section Recombination Parameters:
#' The \code{recpars} parameter contains all parameters required to define the 
#' recombination. \code{recombination_blxAlpha()} understands the following 
#' fields in \code{recpars}:
#' \itemize{
#'    \item \code{alpha} : extrapolation parameter for 'best' parent vector.\cr
#'    Accepts real value \code{0 <= alpha <= 0.5}.
#'    \item \code{beta} : extrapolation parameter for 'worst' parent vector.\cr
#'    Accepts real value \code{0 <= beta <= 0.5}. 
#' }
#' 
#'  @section Warning:
#'  This recombination operator evaluates the candidate solutions in \code{M}, 
#'  which adds an extra \code{popsize} evaluations per iteration.
#'
#' @section References:
#' F. Herrera, M. Lozano, A. M. Sanchez, "A taxonomy for the crossover
#' operator for real-coded genetic algorithms: an experimental study", 
#' International Journal of Intelligent Systems 18(3) 309-338, 2003.
#'
#' @param X population matrix (original)
#' @param M population matrix (mutated) 
#' @param recpars recombination parameters (see \code{Recombination parameters} 
#' for details)
#' 
#' @return Matrix \code{U} containing the recombined population
#' 
#' @export

recombination_blxAlphaBeta <- function(X, M, recpars) {
  
  # Get access to variables in the calling environment
  env <- parent.frame()

  # ========== Error catching and default value definitions
  if (!("alpha" %in% names(recpars))){
    stop("recombination_blxAlphaBeta() requires field alpha in recpars")
  }
  if (!("beta" %in% names(recpars))){
    stop("recombination_blxAlphaBeta() requires field beta in recpars")
  }
  if(!is.numeric(recpars$alpha)){
    stop("recombination_blxAlphaBeta() requires a numeric recpars$alpha")
  }
  if(!is.numeric(recpars$beta)){
    stop("recombination_blxAlphaBeta() requires a numeric recpars$beta")
  }
  if(!(0 <= recpars$alpha & recpars$alpha <= 0.5)){
    stop("recombination_blxAlphaBeta() requires 0 <= recpars$alpha <= 0.5")
  }
  if(!(0 <= recpars$beta & recpars$beta <= 0.5)){
    stop("recombination_blxAlphaBeta() requires 0 <= recpars$beta <= 0.5")
  }
  if (!identical(dim(X), dim(M))) {
    stop("recombination_blxAphaBeta() requires dim(X) == dim(M)")
  }
  if (!all(c("J", "probpars", "nfe") %in% names(env))){
    stop("recombination_blxAphaBeta() requires calling environment to contain 
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
  
  # Get infimum and supremum values, and interval lengths
  Cmin <- pmin(X, M)
  Cmax <- pmax(X, M)
  I    <- Cmax - Cmin
  
  # Get 'best' and 'worst' parents
  C1 <- X * X.is.best + M * !X.is.best
  C2 <- M * X.is.best + X * !X.is.best
  
  S <- (C1 <= C2)
 
  # Return recombined population 
  return(pmin(C1, C2) - 
           I * (recpars$alpha * S + recpars$beta * !S) + 
           randM(X) * (I * ( 1 + recpars$alpha + recpars$beta)))
}
