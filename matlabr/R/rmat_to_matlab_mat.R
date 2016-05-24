#' @title Convert R matrix to matlab matrix
#'
#' @description This function takes in an R matrix then turns it into 
#' a matrix in matlab
#' @param x matrix of values
#' @param matname Object in matlab to be assigned
#' @param transpose Transpose the matrix
#' @export
#' @return Character scalar of matlab code
rmat_to_matlab_mat = function(x, matname = NULL, transpose = FALSE){
  x = as.matrix(x)
  x = apply(x, 1, paste, collapse = ", ")
  x = paste(x, collapse = "; ")
  x = paste0("[", x, "]", ifelse(transpose, "'", ""), ";")
  if (!is.null(matname)) x = paste0(matname, " = ", x)
  x
}