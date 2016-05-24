#' Geometric recombination for DE
#' 
#' Implements the "/geo" (geometric) recombination for the ExpDE framework
#' 
#' @section Recombination Parameters:
#' The \code{recpars} parameter contains all parameters required to define the 
#' recombination. \code{recombination_geo()} understands the following 
#' fields in \code{recpars}:
#' \itemize{
#'    \item \code{alpha} : exponent for geometrical recombination.\cr
#'    Accepts numeric value \code{0 <= alpha <= 1} or \code{NULL} (in which 
#'    case a random value is chosen for each recombination).\cr
#'    Defaults to \code{alpha = 0.5}.
#'}
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

recombination_geo <- function(X, M, recpars = list(alpha = 0.5)) {

  # ========== Error catching and default value definitions
  if (!identical(dim(X), dim(M))) {
    stop("recombination_geo() requires dim(X) == dim(M)")
  }
  if (!("alpha" %in% names(recpars))){
    recpars$alpha <- 0.5
  }
  if (!is.null(recpars$alpha) && !(0 < recpars$alpha && recpars$alpha <= 1)) {
    stop("recombination_geo() requires numeric 0 < recpars$alpha <= 1")
  }
  # ==========
  
  # Get all values to the interval [.25, .75] before performing the recombination
  # (the resulting recombined matrix must be later restored to the original 
  # range)
  mins <- pmin(X, M)
  maxs <- pmax(X, M)
  eps <- 1e-15
  X <- 0.25 + 0.5*(X - mins) / (maxs - mins + eps)
  M <- 0.25 + 0.5*(M - mins) / (maxs - mins + eps)
  
  if(is.null(recpars$alpha)){ # use a random value for each recombination
    alpha <- matrix(rep(stats::runif(nrow(X)),
                        times = ncol(X)),
                    ncol = ncol(X))
  } else{ # use the given (or default) alpha value for all recombinations
    alpha <- recpars$alpha + 0*X
  }
  
  
  # Randomize which parent will use exponent alpha and which will use 
  # (1-alpha)
  inv.alpha <- as.logical(round(stats::runif(nrow(X))))
  alpha[inv.alpha, ] <- 1 - alpha[inv.alpha, ]
  
  # Build recombined population and return it to the original range
  U <- X^alpha * M^(1 - alpha)
  return(mins + (U - 0.25) * (maxs - mins) / 0.5)
}
