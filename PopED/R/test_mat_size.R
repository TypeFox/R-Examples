#' Test to make sure that matricies are the right size
#' 
#' @param correct_size the correct size of a matrix
#' @param mat The matrix to test.
#' @param name The name of the matrix as a string.
#' @example tests/testthat/examples_fcn_doc/examples_test_mat_size.R
#' @export
#' @keywords internal
## Function translated using 'matlab.to.r()'
## Then manually adjusted to make work
## Author: Andrew Hooker

test_mat_size <- function(correct_size,mat,name){
# checks to make sure that matricies are the right size
    
    
    if(all(correct_size==size(mat))){
        return(1)
    } else {
        #stop(sprintf('%s does not have the right dimensions\n Dimensions should be %g X %g', name,correct_size[1],
        #  correct_size[2])) 
        tmp1 <- paste(size(mat),collapse="x")
        tmp2 <- paste(correct_size,collapse="x")
        stop(name, " has dimensions ", tmp1,
             " and should be ", tmp2)
    }
}	

