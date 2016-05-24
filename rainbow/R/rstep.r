rstep <- function (x, FUN = Qn, order = 4, r = matrix.rank(x), mean = TRUE)
{
    if (order < 1)
        stop("Order must be positive.")
    X <- t(x)
    p <- ncol(X)
    n <- nrow(X)
    p1 <- min(order, r, floor(n/2))
    S <- numeric(p1)
    Bnorm <- numeric(n)
    V <- eig <- matrix(0, p, p1)
    Transfo <- diag(p)
    if (mean) {
        med <- L1median2(X, method = "hoss")
        xxx <- xx <- sweep(X, 2, med)
    }
    else xxx <- xx <- X
    for (l in 1:p1) {
        B <- xxx
        for (i in 1:n) Bnorm[i] <- norm(B[i, ], 2)
        Bnormr <- Bnorm[Bnorm > 1e-12]
        B <- B[Bnorm > 1e-12, ]
        A <- diag(1/Bnormr) %*% B
        Y <- xxx %*% t(A)
        s <- colQn(Y)
        j <- order(s, decreasing = TRUE)[1]
        S[l] <- s[j]
        V[l:p, l] <- A[j, ]
        Base <- diag(p - l + 1)
        ndiff <- norm(Base[, 1] - V[l:p, l], Inf)
        if (ndiff > 1e-12) {
            if (sum(V[l:p, l] * Base[, 1]) < 0)
                V[l:p, l] <- -V[l:p, l]
            u <- matrix(Base[, 1] - V[l:p, l], ncol = 1) / c(norm(Base[,
                1] - V[l:p, l]))
            U <- Base - 2 * repmat(t(u) %*% Base, p - l + 1,
                1) * repmat(u, 1, p - l + 1)
        }
        else U <- Base
        eig[, l] <- Transfo %*% V[, l]
        if (l < p1) {
            Edge <- diag(p)
            Edge[l:p, l:p] <- U
            Transfo <- Transfo %*% Edge
            xxx <- xxx %*% U
            xxx <- as.matrix(xxx[, -1])
        }
    }
    coef <- xx %*% eig
    if (mean) {
        basis <- cbind(med, eig)
        coef <- cbind(rep(1, n), coef)
    }
    else basis <- eig
    return(list(basis = basis, coeff = coef, X = xx))
}
