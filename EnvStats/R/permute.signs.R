permute.signs <-
function (n) 
{
    if (!is.vector(n, mode = "numeric") || is.factor(n) || length(n) != 
        1 || n != trunc(n) || n <= 0) 
        stop("'n' must be a positive integer")
    vec <- c(-1, 1)
    nr <- 2^n
    mat <- matrix(NA, ncol = n, nrow = nr)
    for (i in 1:n) {
        x <- rep(vec, each = 2^(i - 1))
        mat[, i] <- rep(x, length = nr)
    }
    mat
}
