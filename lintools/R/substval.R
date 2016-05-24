#' Substitute a value in a system of linear (in)equations
#'
#' @section Details:
#' A system of the form \code{Ax <= b} can be simplified if one or more of the
#' \code{x[i]} values is fixed.
#'
#'
#' @param A \code{[numeric]} matrix
#' @param b \code{[numeric]} vector 
#' @param variables \code{[numeric|logical|character]} vector of column indices in \code{A}
#' @param values \code{[numeric]} vecor of values to substitute.
#' @param remove_columns \code{[logical]} Remove spurious columns when substituting?
#'
#' @return A \code{list} with the following components:
#' 
#' \itemize{
#'   \item{\code{A}: the \code{A} corresponding to the simplified sytem.}
#'   \item{\code{b}: the constant vector corresponding to the new system}
#' }
#'
#'
#' @export
subst_value <- function(A, b, variables, values, remove_columns=FALSE){
  check_sys(A=A, b=b)
  if ( is.character(variables) ){
    variables <- match(variables,colnames(A))
  }

  b <- as.vector(b - A[,variables,drop=FALSE] %*% values)
  if (remove_columns){
    A <- A[,-variables,drop=FALSE]
  } else {
    A[,variables] <- 0
  }
  list(A=A, b=b)
}
