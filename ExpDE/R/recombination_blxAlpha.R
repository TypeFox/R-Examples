#' Blend Alpha recombination for DE
#' 
#' Implements the "/blxAlpha" (Blend Alpha) recombination for the ExpDE 
#' framework
#' 
#' @section Recombination Parameters:
#' The \code{recpars} parameter contains all parameters required to define the 
#' recombination. \code{recombination_blxAlpha()} understands the following 
#' fields in \code{recpars}:
#' \itemize{
#'    \item \code{alpha} : extrapolation parameter.\cr
#'    Accepts real value \code{0 <= alpha <= 0.5}.
#' }
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

recombination_blxAlpha <- function(X, M, recpars) {

  # ========== Error catching and default value definitions
  if (!("alpha" %in% names(recpars))){
    stop("recombination_blxAlpha() requires field alpha in recpars")
  }
  if(!is.numeric(recpars$alpha)){
    stop("recombination_blxAlpha() requires a numeric rectpars$alpha")
  }
  if(!(0 <= recpars$alpha & recpars$alpha <= 0.5)){
    stop("recombination_blxAlpha() requires 0 <= recpars$alpha <= 0.5")
  }
  if (!identical(dim(X), dim(M))) {
    stop("recombination_blxApha() requires dim(X) == dim(M)")
  }
  # ==========
  
  Cmin <- pmin(X, M)
  Cmax <- pmax(X, M)
  I    <- Cmax - Cmin
 
  # Return recombined population
  return (Cmin - I * recpars$alpha + 
            randM(X) * I * (1 + 2 * recpars$alpha))
}
