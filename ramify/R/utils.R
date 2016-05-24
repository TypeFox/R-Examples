#' Shorten a Vector
#' 
#' Shorten a vector using \code{...} notation.
#' 
#' @keywords internal
add_dots <- function(x, pos = 3) {
  if (length(x) >= pos + 2) {
    c(x[seq_len(pos-1)], "...", x[length(x)])
  } else {
    x
  }
}


#' Describe a Matrix
#' 
#' Prints a short description about a matrix.
#' 
#' @keywords internal
desc_mat <- function(x) {
  # paste(paste(dim(x), collapse = " by "), "matrix of", paste0(typeof(x), "s"))
  paste(paste(dim(x), collapse = " x "), "matrix of", paste0(typeof(x), "s:"))
}
