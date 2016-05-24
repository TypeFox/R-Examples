#' Test if any matrix element is above a max value. 
#' 
#' Function tests is any matrix element is above a maximum value.  For those elements 
#' the function sets those values to the maximum value.
#'  
#' 
#' @param mat A matrix.
#' @param max_mat A matrix the same size as mat with the  maximum allowed value of that element.
#' @return A matrix  
#' @family matrix_manipulation
#' @example tests/testthat/examples_fcn_doc/examples_test_for_max.R
#' @export
#' @keywords internal
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

test_for_max <- function(mat,max_mat){

if(any(any(mat>max_mat))){
  for(i in 1:size(mat,1)){
    for(j in 1:size(mat,2)){
      if(mat[i,j]>max_mat[i,j]){
        mat[i,j]=max_mat[i,j]
      }
    }
  }
}
ret=mat
return( ret)
}
