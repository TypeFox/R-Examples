##
##  t r i q u a d . R  Gaussian Triangle Quadrature
##


triquad <- function(f, x, y, n = 10, tol = 1e-10, ...) {
    if (!is.numeric(x) || length(x) != 3 ||
        !is.numeric(y) || length(y) != 3)
        stop("Arguments 'x' and 'y' must be numeric vectors of length 3.")
    n <- floor(abs(n))
    if (n <= 2)
        stop("Argument 'n' must be an integer greater or equal 2.")

    fun <- match.fun(f)
    f <- function(x, y) fun(x, y, ...)

    v <- cbind(x, y)
    rel.tol <- Inf
    kmax <- 5
    k <- 1
    while (rel.tol > tol && k <= kmax) {
        G1 <- .tricoef(v, N = n)
        I1 <- t(G1$Wx) %*% f(G1$X, G1$Y) %*% G1$Wy
        G2 <- .tricoef(v, N = 2*n+1)
        I2 <- t(G2$Wx) %*% f(G2$X, G2$Y) %*% G2$Wy
        rel.tol <- abs(I1 - I2)
        k <- k + 1
    }
    return(c(I2))
}


.tricoef <- function(v, N = 32, ...) {
    eps <- .Machine$double.eps

    n <- 1:N
    nnk <- 2*n + 1
    A <- c(1/3, repmat(1,1,N) / (nnk * (nnk+2)))
    
    n <- 2:N
    nnk <- nnk[n]
    B1 <- 2/9
    nk <- n+1
    nnk2 <- nnk * nnk
    B <- 4*(n * nk)^2 / (nnk2 * nnk2 - nnk2)
    ab <- cbind(A, c(2, B1, B))
    s <- sqrt(ab[2:N, 2])

    VX <- eigen(Diag(ab[1:N,1], 0) + Diag(s, -1) + Diag(s, 1))
    X <- VX$values
    V <- VX$vectors
    I <- order(X)
    X <- X[I]
    x <- (X + 1) / 2
    wx <- ab[1,2] * as.matrix(V[1,I])^2 / 4

    N <- N-1; N1 <- N+1; N2 <- N+2
    y <- cos((2*as.matrix(N:0)+1)*pi/(2*N+2))

    L <- zeros(N1, N2)
    y0=2
    iter=0;
    while (max(abs(y-y0)) > eps) {
        L[, 1] <- 1
        L[, 2] <- y  
        for (k in 2:N1) {
            L[, k+1] <- ( (2*k-1) * y * L[, k] - (k-1) * L[,k-1] ) / k
        }
        Lp <- N2 * ( L[, N1] - y * L[,N2] ) / (1-y^2) 
        y0 <- y
        y <- y0 - L[, N2]/Lp
        iter <- iter+1
    }
    cc <- matrix(c(1,0,0, -1,0,1, 0,1,-1), 3, 3, byrow = TRUE) %*% v
    t1 <- (1+y)/2
    Wx <- abs(det(cc[2:3,])) * wx
    Wy <- 1/((1-y^2) * Lp^2) * (N2/N1)^2
    mg <- meshgrid(t1, x)
    t2 <- mg$X; xx <- mg$Y 
    yy <- t2 * xx
    X  <- cc[1,1] + cc[2,1]*xx + cc[3,1]*yy
    Y  <- cc[1,2] + cc[2,2]*xx + cc[3,2]*yy

    return(list(X = X, Y = Y, Wx = Wx, Wy = Wy))
}
