#' /rand mutation for DE
#' 
#' Implements the "/rand/nvecs" mutation for the ExpDE framework
#' 
#' @section Mutation Parameters:
#' The \code{mutpars} parameter contains all parameters required to define the 
#' mutation. \code{mutation_rand()} understands the following fields in 
#' \code{mutpars}:
#' \itemize{
#'    \item \code{f} : scaling factor for difference vector(s).\cr
#'    Accepts numeric vectors of size 1 or \code{nvecs}.
#'    \item \code{nvecs} : number of difference vectors to use.\cr 
#'        Accepts \code{1 <= nvecs <= (nrow(X)/2 - 2)}\cr
#'        Defaults to 1.
#' }
#' 
#' @param X population matrix
#' @param mutpars mutation parameters (see \code{Mutation parameters} for details)
#' 
#' @return Matrix \code{M} containing the mutated population
#' 
#' @section References:
#' K. Price, R.M. Storn, J.A. Lampinen, "Differential Evolution: A 
#' Practical Approach to Global Optimization", Springer 2005
#' 
#' @export

mutation_rand <- function(X, mutpars){

  # ========== Error catching and default value definitions
  if (!("nvecs" %in% names(mutpars))) mutpars$nvecs <- 1
  if (!(mutpars$nvecs %in% 1:(nrow(X)/2 - 2))){
    stop("mutation_rand() requires integer 1 <= mutpar$nvecs <= (popsize/2 - 2)")
  }
  if (!("f" %in% names(mutpars))){
    stop("mutation_rand() requires field f in mutpars")
  }
  if (length(mutpars$f) == 1) mutpars$f <- rep(mutpars$f, 
                                               mutpars$nvecs)
  # ==========
  
  # Matrix indices for mutation (r1 != r2 != r3 != ... != rn)
  R <- lapply(X = rep(nrow(X), 
                      times = nrow(X)),
              FUN = sample.int,
              size = 1 + 2*mutpars$nvecs,
              replace = FALSE)

    
  # Auxiliary function: make a single mutation
  randmut <- function(pos, Pop, f){
    Xr1   <- Pop[pos[1], ]
    diffs <- matrix(pos[-1],
                    ncol=2,
                    byrow=TRUE)
    if (nrow(diffs) == 1) {
      wdiffsum <- f*(Pop[diffs[, 1], ] - Pop[diffs[, 2], ])
    } else {
      wdiffsum <- colSums(f * (Pop[diffs[, 1], ] - Pop[diffs[, 2], ]))
    }
    return(Xr1 + wdiffsum)
  }
  
  # Apply mutation
  M <- lapply(R, 
              FUN = randmut, 
              Pop = X, 
              f = mutpars$f)
  
  return(matrix(data  = unlist(M), 
                nrow  = nrow(X), 
                byrow = T))
}