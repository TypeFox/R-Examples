#' Exponential recombination for DE
#' 
#' Implements the "/exp" (exponential) recombination for the ExpDE framework
#' 
#' @section Recombination Parameters:
#' The \code{recpars} parameter contains all parameters required to define the 
#' recombination. \code{recombination_exp()} understands the following 
#' fields in \code{recpars}:
#' \itemize{
#'    \item \code{cr} : component-wise probability of selection as a cut-point.
#'    \cr
#'    Accepts numeric value \code{0 < cr <= 1}.
#' }
#' 
#' @section References:
#' K. Price, R.M. Storn, J.A. Lampinen, "Differential Evolution: A 
#' Practical Approach to Global Optimization", Springer 2005
#'
#' @param X population matrix (original)
#' @param M population matrix (mutated) 
#' @param recpars recombination parameters (see \code{Recombination parameters} 
#' for details)
#' 
#' @return Matrix \code{U} containing the recombined population
#' 
#' @export

recombination_exp <- function(X, M, recpars) {
  
  # ========== Error catching and default value definitions
  if (!identical(dim(X), dim(M))) {
    stop("recombination_exp() requires dim(X) == dim(M)")
  }
  if (!("cr" %in% names(recpars))){
    stop("recombination_exp() requires field cr in recpars")
  }
  if (!(0 < recpars$cr & recpars$cr <= 1)) {
    stop("recombination_exp() requires numeric 0 < recpars$cr <= 1")
  }
  # ==========
  
  # Start points for mutation: for each row, a value between 1 and length(x),
  # uniformly distributed
  mut.start <- sample.int(n       = ncol(X),
                          size    = nrow(X),
                          replace = TRUE)
  
  # End points for mutation: for each row, a value between mut.start and 
  # (mut.start + length(x) - 1), exponentially distributed
  probs <- recpars$cr^(1:ncol(X) - 1) - recpars$cr^(1:ncol(X))
  mut.end    <- mut.start + sample(x    = 1:ncol(X) - 1,
                                   size = nrow(X),
                                   replace = TRUE,
                                   prob = probs / sum(probs))
  
  # Helper function for setting mutation indices: 
  # for each row wrap around the end of the vector, 
  # e.g., if n = 5, s = 3 and e = 6, returns z = [1, 0, 1, 1, 1] (pos 3,4,5,1)
  # e.g., if n = 5, s = 1 and e = 1, returns z = [1, 0, 0, 0, 0] (pos 1)
  setfun <- function(n, s, e) {
    z<-numeric(n)
    z[(s - 1):(e - 1) %% n + 1] <- 1
    z
  }
  
  # Recombination matrix - using mapply() to apply over multiple indexed objects
  R <- t(mapply(FUN      = setfun, 
                n        = rep(ncol(X), nrow(X)),
                s        = mut.start,
                e        = mut.end,
                SIMPLIFY = TRUE))
  
  # Return recombined population
  return(R*M + (1 - R)*X)
}
