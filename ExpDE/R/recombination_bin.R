#' /bin recombination for DE
#' 
#' Implements the "/bin" (binomial) recombination for the ExpDE framework
#' 
#' @section Recombination Parameters:
#' The \code{recpars} parameter contains all parameters required to define the 
#' recombination. \code{recombination_bin()} understands the following fields in 
#' \code{recpars}:
#' \itemize{
#'    \item \code{cr} : component-wise probability of using the value in 
#'                      \code{M}.\cr
#'                      Accepts numeric value \code{0 < cr <= 1}.
#'    \item \code{minchange} : logical flag to force each new candidate solution 
#'                             to inherit at least one component from its 
#'                             mutated 'parent'.\cr
#'                             Defaults to TRUE
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

recombination_bin <- function(X, M, recpars) {

  # ========== Error catching and default value definitions
  if (!("cr" %in% names(recpars))){
    stop("recombination_bin() requires field cr in recpars")
  }
  if (!(0 < recpars$cr & recpars$cr <= 1)) {
    stop("recombination_bin() requires numeric 0 < recpars$cr <= 1")
  }
  if (!identical(dim(X),dim(M))) {
    stop("recombination_bin() requires dim(X) == dim(M)")
  }
  if (!("minchange" %in% names(recpars))){
    recpars$minchange <- TRUE
  }
  # ==========
  
  # Recombination matrix
  R <- randM(X) < recpars$cr
  
  if (recpars$minchange){
    indx    <- which(rowSums(R) == 0)
    cor.mat <- cbind(indx, 
                     sample.int(n       = ncol(X),
                                size    = length(indx), 
                                replace = TRUE))
    R[cor.mat[,1],cor.mat[,2]] <- TRUE
  }
  
  # Return recombined population
  return(R*M + (1 - R)*X)
}
