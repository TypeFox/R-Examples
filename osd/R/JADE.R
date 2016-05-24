# JADE <-
# function (X, n.comp = NULL, eps = 1e-06, maxiter = 100, na.action = na.fail, center=T) 
# {
    # X <- na.action(X)
    # if (!all(sapply(X, is.numeric))) 
        # stop("'X' must be numeric")
    # X <- as.matrix(X)
    # data.matrix <- X
    # N <- dim(X)[1]
    # X.cols <- dim(X)[2]
    # if (is.null(n.comp)) 
        # n.comp <- X.cols
    # if (n.comp > X.cols) 
        # stop("'n.comp' must be smaller than number of columns of X")
    # if(center==T) X <- scale(X, scale = F)
    # Col.center <- attr(X, "scaled:center")
    # data.X <- X
    # eigen.X <- eigen(crossprod(X, X)/N)
    # U1 <- eigen.X$vectors
    # D1 <- eigen.X$values
    # puiss <- sort(D1, decreasing = F)
    # k <- 1:length(D1)
    # U <- U1[, rev(k)]
    # rangeW <- (X.cols - n.comp + 1):X.cols
    # scales <- sqrt(puiss[rangeW])
    # if (n.comp > 1) {
        # W <- diag(1/scales) %*% t(U[, rangeW])
        # iW <- U[, rangeW] %*% diag(scales)
    # }
    # else {
        # W <- t(U[, rangeW])/scales
        # iW <- U[, rangeW] * scales
    # }
    # X <- X %*% t(W)
    # dimsymm <- (n.comp * (n.comp + 1))/2
    # nbcm <- dimsymm
    # CM <- matrix(0, nrow = n.comp * nbcm, ncol = n.comp)
    # R <- diag(n.comp)
    # Qij <- matrix(0, ncol = n.comp, nrow = n.comp)
    # Xim <- numeric(n.comp)
    # Xjm <- numeric(n.comp)
    # scale2 <- rep(1, n.comp)/N
    # Range <- 1:n.comp
    # for (im in 1:n.comp) {
        # Xim <- X[, im]
        # Qij <- t((Xim * Xim) %*% t(scale2) * X) %*% X - R - 2 * 
            # R[, im] %*% t(R[, im])
        # CM[Range, ] = Qij
        # Range <- Range + n.comp
        # if (im > 1) {
            # for (jm in (1:(im - 1))) {
                # Xjm <- X[, jm]
                # Qij <- t((Xim * Xjm) %*% t(scale2) * X) %*% X - 
                  # R[, im] %*% t(R[, jm]) - R[, jm] %*% t(R[, 
                  # im])
                # CM[Range, ] = sqrt(2) * Qij
                # Range <- Range + n.comp
            # }
        # }
    # }
    # V <- t(rjd.fortran(CM)$V)
    # B <- V %*% W
    # S <- tcrossprod(data.X, B)
    # kurt <- rep(0, n.comp)
    # for (j in 1:n.comp) {
        # kurt[j] <- mean(S[, j]^4) - 3
    # }
    # P <- matrix(0, n.comp, n.comp)
    # ord <- order(kurt, decreasing = TRUE)
    # for (j in 1:n.comp) {
        # P[j, ord[j]] <- 1
    # }
    # B <- P %*% B
    # if (n.comp > 1) {
        # B <- diag(sign(rowMeans(B))) %*% B
    # }
    # else B <- sign(sum(B)) * B
    # if (n.comp == X.cols) 
        # A <- solve(B)
    # else A <- t(B) %*% solve(B %*% t(B))
    # S <- tcrossprod(data.X, B)
    # colnames(S) <- paste("IC.", 1:n.comp, sep = "")
    # res <- list(A = A, W = B, S = S, Xmu = Col.center)
    # class(res) <- "bss"
    # return(res)
# }

JADE <- function (X, n.comp = NULL, eps = 1e-06, maxiter = 100, na.action = na.fail, center=T) 
{
    X <- na.action(X)
    if (!all(sapply(X, is.numeric))) 
        stop("'X' must be numeric")
    X <- as.matrix(X)
    data.matrix <- X
    N <- dim(X)[1]
    X.cols <- dim(X)[2]
    if (is.null(n.comp)) 
        n.comp <- X.cols
    if (n.comp > X.cols) 
        stop("'n.comp' must be smaller than number of columns of X")
    if(center==T) X <- scale(X, scale = F)
    Col.center <- attr(X, "scaled:center")
    data.X <- X
    eigen.X <- eigen(crossprod(X)/N, symmetric = TRUE)
    U1 <- eigen.X$vectors
    D1 <- eigen.X$values
    puiss <- sort(D1, decreasing = FALSE)
    k <- 1:length(D1)
    U <- U1[, rev(k)]
    rangeW <- (X.cols - n.comp + 1):X.cols
    scales <- sqrt(puiss[rangeW])
    if (n.comp > 1) {
        W <- tcrossprod(diag(1/scales), U[, rangeW])
        iW <- tcrossprod(U[, rangeW], diag(scales))
    }
    else {
        W <- t(U[, rangeW])/scales
        iW <- U[, rangeW] * scales
    }
    X <- tcrossprod(X, W)
    dimsymm <- (n.comp * (n.comp + 1))/2
    nbcm <- dimsymm
    CM <- matrix(0, nrow = n.comp * nbcm, ncol = n.comp)
    R <- diag(n.comp)
    Qij <- matrix(0, ncol = n.comp, nrow = n.comp)
    Xim <- numeric(n.comp)
    Xjm <- numeric(n.comp)
    scale2 <- rep(1, n.comp)/N
    Range <- 1:n.comp
    for (im in 1:n.comp) {
        Xim <- X[, im]
        Qij <- crossprod((tcrossprod((Xim * Xim), scale2) * X), 
            X) - R - 2 * tcrossprod(R[, im])
        CM[Range, ] <- Qij
        Range <- Range + n.comp
        if (im > 1) {
            for (jm in (1:(im - 1))) {
                Xjm <- X[, jm]
                Qij <- crossprod((tcrossprod((Xim * Xjm), scale2) * 
                  X), X) - tcrossprod(R[, im], R[, jm]) - tcrossprod(R[, 
                  jm], R[, im])
                CM[Range, ] <- sqrt(2) * Qij
                Range <- Range + n.comp
            }
        }
    }
    V <- t(frjd(CM)$V)
    B <- V %*% W
    S <- tcrossprod(data.X, B)
    kurt <- rep(0, n.comp)
    for (j in 1:n.comp) {
        kurt[j] <- mean(S[, j]^4) - 3
    }
    P <- matrix(0, n.comp, n.comp)
    ord <- order(kurt, decreasing = TRUE)
    for (j in 1:n.comp) {
        P[j, ord[j]] <- 1
    }
    B <- P %*% B
    if (n.comp > 1) {
        B <- sweep(B, 1, sign(rowMeans(B)), "*")
    }
    else B <- sign(sum(B)) * B
    if (n.comp == X.cols) 
        A <- solve(B)
    else A <- crossprod(B, solve(tcrossprod(B)))
    S <- tcrossprod(data.X, B)
    colnames(S) <- paste("IC.", 1:n.comp, sep = "")
    res <- list(A = A, W = B, S = S, Xmu = Col.center)
    class(res) <- "bss"
    return(res)
}






