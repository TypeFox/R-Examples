#' Min Max Arithmetical recombination for DE
#' 
#' Implements the "/mmax" (min-max-arithmetical) recombination for the ExpDE 
#' framework
#'
#' @section Warning:
#' This recombination operator evaluates \code{4*popsize} candidate solutions 
#' per iteration of the algorithm. The value of the \code{nfe} counter and the 
#' vector of performance values \code{G} are updated in the calling environment.
#' 
#' @section Recombination Parameters:
#' The \code{recpars} parameter contains all parameters required to define the 
#' recombination. \code{recombination_pbest()} understands the following 
#' fields in \code{recpars}:
#' \itemize{
#'    \item \code{lambda} : Recombination multiplier.\cr
#'                          Optional. Defaults to \code{NULL}
#'                          Accepts numeric value \code{0 < lambda < 1} or 
#'                          \code{NULL} (in which case a random value is 
#'                          independently used for each variable of each 
#'                          recombination pair).
#'}
#' 
#' @section References:
#' F. Herrera, M. Lozano, A. M. Sanchez, "A taxonomy for the crossover
#' operator for real-coded genetic algorithms: an experimental study", 
#' International Journal of Intelligent Systems 18(3):309-338, 2003.\cr
#' F Herrera, M. Lozano,  J.L. Verdegay, "Tuning fuzzy logic controllers by 
#' genetic algorithms.", International Journal of Approximate Reasoning 
#' 12(3):299-315, 1995. \cr
#'
#' @param X population matrix (original)
#' @param M population matrix (mutated) 
#' @param recpars recombination parameters (see \code{Recombination parameters} 
#' for details)
#' 
#' @return Matrix \code{U} containing the recombined population
#' 
#' @export

recombination_mmax <- function(X, M, recpars = list(lambda = NULL)) {
  
  # Get access to variables in the calling environment
  env <- parent.frame()
  
  # ========== Error catching and default value definitions
  if (!identical(dim(X), dim(M))) {
    stop("recombination_mmax() requires dim(X) == dim(M)")
  }
  if (!all(c("probpars", "nfe") %in% names(env))){
    stop("recombination_mmax() requires calling environment to contain 
         variables nfe and probpars")
  }
  if ("lambda" %in% names(recpars)){
    if(!is.null(recpars$lambda)){
      if(!is.numeric(recpars$lambda) || 
         !(0 < recpars$lambda & recpars$lambda < 1) ||
         length(recpars$lambda) != 1){
        stop("recombination_mmax() requires recpars$lambda to be either NULL 
                or a single value between 0 and 1")
      }
    }
  } else recpars$lambda <- NULL
  # ==========
  
  # Define lambda factors
  if(is.null(recpars$lambda)) {
    lambda <- randM(X)
  } else {
    lambda <- matrix(recpars$lambda, 
                     nrow = nrow(X),
                     ncol = ncol(X))
  } 
  
  # Generate trial offspring
  H1 <- lambda * X + (1 - lambda) * M
  H2 <- (1 - lambda) * X + lambda * M
  H3 <- pmin(X, M)
  H4 <- pmax(X, M)
  
  # Evaluate trial offspring
  f1 <- evaluate_population(probpars = env$probpars, 
                            Pop      = H1)
  
  f2 <- evaluate_population(probpars = env$probpars, 
                            Pop      = H2)
  
  f3 <- evaluate_population(probpars = env$probpars, 
                            Pop      = H3)
  
  f4 <- evaluate_population(probpars = env$probpars, 
                            Pop      = H4)
  
  env$nfe <- env$nfe + 4 * nrow(X)
  
  # Get 'winning' offspring
  fbest <- pmin(f1, f2, f3, f4)
  
  # Update performance vector in calling environment
  env$G[f1 == fbest] <- f1[f1 == fbest]
  env$G[f2 == fbest] <- f2[f2 == fbest]
  env$G[f3 == fbest] <- f3[f3 == fbest]
  env$G[f4 == fbest] <- f4[f4 == fbest]
  
  # Assemble output population
  Pop.trialx <- NA * X
  Pop.trialx[f1==fbest, ] <- H1[f1==fbest, ]
  Pop.trialx[f2==fbest, ] <- H2[f2==fbest, ]
  Pop.trialx[f3==fbest, ] <- H3[f3==fbest, ]
  Pop.trialx[f4==fbest, ] <- H4[f4==fbest, ]
  
  # Return population
  return (Pop.trialx)
  }
