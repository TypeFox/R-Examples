#' @title Relabels a vector to consecutive labels
#'
#' @description 
#' This function relabels a vector to have consecutive - no missing in between -
#' labels.  Labels always start at \eqn{1} and increase by one.
#' 
#' For example, \code{c(2, 2, 5)} gets relabeled to \code{c(1, 1, 2)}.
#' 
#' @param vec vector with labels
#' @param order logical; if \code{TRUE} then new state labels are assigned by
#' decreasing number of points in that state. That is, state ``1'' has the most
#' points in the state, followed by state ``2'' etc. 
#' @keywords manip array list arith
#' @export
#' @examples
#' 
#' TempVec = c(10,2,1,2,2,2,10)
#' print(relabel_vector(TempVec))
#' 
#' print(relabel_vector(c(2, 2, 5)))
#' 

relabel_vector <- function(vec, order = FALSE) {
  
  vec_unique <- unique(vec)
  num.elements <- length(vec_unique)
  
  if (order) {
    vec_unique <- vec_unique[order(table(vec), decreasing = TRUE)]
  }
  # TODO: improve this with a vectorized approach
  vec_new <- rep(NA, length(vec))
  for (ii in seq_len(num.elements)) {
    vec_new[vec == vec_unique[ii]] <- ii
  }
  
  invisible(vec_new)
} 
