#' Test if any matrix element is above a minimum value. 
#' 
#' Function tests is any matrix element is above a minimum value.  For those elements 
#' the function sets those values to the minimum value.
#'  
#' 
#' @param mat A matrix.
#' @param min_mat A matrix the same size as mat with the  minimum allowed value of that element.
#' @return A matrix 
#' @family matrix_manipulation 
#' @example tests/testthat/examples_fcn_doc/examples_test_for_min.R
#' @export
#' @keywords internal
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

test_for_min <- function(mat,min_mat){

if(any(any(mat<min_mat))){
  for(i in 1:size(mat,1)){
    for(j in 1:size(mat,2)){
      if(mat[i,j]<min_mat[i,j]){
        mat[i,j]=min_mat[i,j]
      }
    }
  }
}
ret=mat
return( ret)
}
