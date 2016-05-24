#' Head and tail in two dimensions 
#' 
#' Restrict not only the number of rows, but also the number of columns.
#' @param x data.frame or matrix
#' @param m integer. Number of rows to keep.
#' @param n integer. Number of columns to keep.
#' 
#' @author kaiyin
#' @export
head2 = function(x, m=6, n=NULL) {
    if(is.null(n)) {
        n = m
    }
    x[1:m, 1:n]
}

#' @rdname head2
#' @export 
tail2 = function(x, m = 6, n = NULL) {
    if(is.null(n)) {
        n = m
    }
    e1 = nrow(x)
    e2 = ncol(x)
    i1 = (e1 - m + 1)
    i2 = (e2 - n + 1)
    x[i1:e1, i2:e2]
}
