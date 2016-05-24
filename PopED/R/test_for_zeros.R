#' Test if any matrix element is zero. 
#' 
#' Function tests is any matrix element is zero.  For those elements 
#' the function sets those values to the minimum value allowed (not zero).
#' This is used to avoid numerical problems in the FIM calculation.
#'  
#' 
#' @param mat A matrix.
#' @param ourzero A matrix the same size as mat with the  value that zero should be reassigned to.
#' @return A matrix 
#' @family matrix_manipulation 
#' @example tests/testthat/examples_fcn_doc/examples_test_for_zeros.R
#' @export
#' @keywords internal
## Function translated automatically using 'matlab.to.r()'
## Author: Andrew Hooker

test_for_zeros <- function(mat,ourzero){

if(any(any(mat==0))){
  for(i in 1:size(mat,1)){
    for(j in 1:size(mat,2)){
      if(mat[i,j]==0){
        mat[i,j]=ourzero
      }
    }
  }
}
ret=mat
return( ret)
}
