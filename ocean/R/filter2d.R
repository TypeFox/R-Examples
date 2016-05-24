filter.vec <- function(vec, r, i) {
    vec.out <- vec
    for(i in seq(length(vec))) {
        filt <- dnorm(seq(length(vec)), i, r)
        filt <- filt / sum(filt)
        vec.out[i] <- sum(vec * filt)
    }
    return(vec.out)
}

filter.row <- function(mat, r, i) {
    mat[i,] <- filter.vec(mat[i,], r)
    return(mat)
}

filter.col <- function(mat, r, i)
    return(t(filter.row(t(mat), r, i)))

#' Applies a 2d Gaussian filter to \code{mat} with a standard deviation of r
#' cells.
#'
#' Applies a Gaussian blur to a 2D matrix. The matrix is first convoluted
#' with the filter along rows, then along columns.
#'
#' @param mat The matrix to filter.
#' @param r The standard deviation of the filter in matrix cells.
filter2d <- function(mat, r) {
    for(i in seq(nrow(mat)))
        mat <- filter.row(mat, r, i)
    for(i in seq(ncol(mat)))
        mat <- filter.col(mat, r, i)
    return(mat)
}
