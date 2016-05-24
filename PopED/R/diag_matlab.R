#' Function written to match MATLAB's diag function
#' 
#' There are some differences between tha MATLAB and the R version of diag.
#' Specifically, if a 1xN or a Nx1 matrix is supplied to the R
#' \code{\link{diag}} function then just the first element of this vector is
#' returned. This function tries to match the MATLAB version in handling vectors
#' (matricies with one dimension equal to one), and will return a diagonal
#' matrix in these situations.
#' 
#' @param mat Either a vector to make into a diagonal matrix or a matrix you 
#'   want to extract the diagonal from
#' @return Either a diagonal matrix or the diagonal of a matrix.
#' @family MATLAB
#' @family matrix_manipulation
#' @example tests/testthat/examples_fcn_doc/examples_diag_matlab.R
#' @export
#' @keywords internal

diag_matlab <- function(mat){
    dim.mat <- dim(mat)
    if(!is.null(dim.mat)){
        if(any(dim.mat==1)){
            if(!all(dim.mat==1)){                
                mat <- mat[,,drop=T]
            }
        }
    }
    return(diag(mat))
}
