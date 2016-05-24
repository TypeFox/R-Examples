repmat = function (A, m, n = m)
{
    A <- as.matrix(A)
    tmp <- matrix(rep(t(A), m), nrow = m * nrow(A), byrow = TRUE)
    return(matrix(rep(tmp, n), ncol = n * ncol(tmp)))
}
