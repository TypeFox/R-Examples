circshift <- function(a, sz) {
    if (is.null(a)) return(a)
    
    if (is.vector(a) && length(sz) == 1) {
        n <- length(a)
        s <- sz %% n
        a <- a[(1:n-s-1) %% n + 1]

    } else if (is.matrix(a) && length(sz) == 2) {
        n <- nrow(a); m <- ncol(a)
        s1 <- sz[1] %% n
        s2 <- sz[2] %% m
        a <- a[(1:n-s1-1) %% n + 1, (1:m-s2-1) %% m + 1]
    } else
        stop("Length of 'sz' must be equal to the no. of dimensions of 'a'.")

    return(a)
}
