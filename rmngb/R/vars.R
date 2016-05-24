colVars <- function(x, na.rm = FALSE, dims = 1, unbiased = TRUE,
                    SumSquares = FALSE, twopass = TRUE) {
    if (SumSquares) return(colSums(x^2, na.rm, dims))
    
    N <- colSums(!is.na(x), FALSE, dims)
    Nm1 <- if (unbiased) N - 1 else N
    
    if (twopass) {
        x <- if (dims == length(dim(x))) x - mean(x, na.rm = na.rm) else
            sweep(x, (dims + 1):length(dim(x)), colMeans(x, na.rm, dims))
    }
    (colSums(x^2, na.rm, dims) - colSums(x, na.rm, dims)^2 / N) / Nm1
}

rowVars <- function(x, na.rm = FALSE, dims = 1, unbiased = TRUE,
                    SumSquares = FALSE, twopass = TRUE) {
    if (SumSquares) return(rowSums(x^2, na.rm, dims))
    
    N <- rowSums(!is.na(x), FALSE, dims)
    Nm1 <- if (unbiased) N - 1 else N
    
    if (twopass) {
        x <- if (dims == 0) x - mean(x, na.rm = na.rm) else
            sweep(x, 1:dims, rowMeans(x, na.rm, dims))
    }
    (rowSums(x^2, na.rm, dims) - rowSums(x, na.rm, dims)^2 / N) / Nm1
}
