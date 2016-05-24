##
##  e x p m . R  Matrix Exponential
##


expm <- function(A, np = 128) {
    if (!is.numeric(A) || !is.matrix(A) || nrow(A) != ncol(A))
        stop("Argument 'A' must be a square numeric matrix.")
    if (!is.numeric(np) || length(np) != 1 ||
        floor(np) != ceiling(np) || np < 2)
        stop("Argument 'np' must be an integer greater or equal to 2.")

    N <- nrow(A)
    circle <- exp(2i*pi*(1:np)/np)      # generate np unit roots
    z0 <- ceiling(mean(range(Re(eig(A)))) + 0.1)
    radius <- ceiling(max(abs(eig(A) - z0)) + 0.1)
    z <- z0 + radius*circle

    I <- eye(N); B <- zeros(N)
    for (i in 1:np) {
      R <- inv(z[i]*I - A)                # resolvent matrix at point z(i)
      B <- B + R * (z[i]-z0) * exp(z[i])  # add up contributions to integral
    }

    B <- Re(B)/np
    return(zapsmall(B))
}


logm <- function(A) {
    if (!is.numeric(A) || !is.matrix(A) || nrow(A) != ncol(A))
        stop("Argument 'A' must be a square numeric matrix.")

    E <- eigen(A)
    e <- E$values
    if (any(Im(e) == 0 && Re(e) <= 0))
        stop("A must not have any nonpositive real eigenvalues.")


    D <- diag(log(E$values))
    X <- E$vectors %*% D %*% solve(E$vectors)
    return(Re(X))
}