#' Reshape a vector into a matrix
#' 
#' This function reshapes a vector into a matrix whose row elements are taken from the vector. Orders of elements keep unchanged
#' from the vector.
#' 
#' @param x a vector or a time series
#' @param w a number specifying number of columns of the matrix
#' @return a matrix
#' @export
#' @examples
#' x=ts2mat(c(1:(128*20)),128)
#' dim(x)
#' x[1,1:20]


ts2mat <- function(x, w) {
    m = length(x)
    mat <- matrix(0, floor(m/w), w)
    for (i in 1:floor(m/w)) {
        mat[i, ] = x[((i - 1) * w + 1):(i * w)]
    }
    return(mat)
    
} 
