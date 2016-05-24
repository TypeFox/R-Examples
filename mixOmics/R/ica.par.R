# The function ica.par and ica.def are borrowed from the fastICA package (see references).

ica.par <- function (X, ncomp, tol, fun, alpha, max.iter, verbose, w.init)
{
    Diag <- function(d) if(length(d) > 1L) diag(d) else as.matrix(d)
    n <- nrow(X)
    p <- ncol(X)
    W <- w.init+matrix(c(1:9),ncomp,ncomp)
    sW <- La.svd(W)
    W <- sW$u %*% Diag(1/sW$d) %*% t(sW$u) %*% W
    W1 <- W
    lim <- rep(1000, max.iter)
    it <- 1
    if (fun == "logcosh") {
        while (lim[it] > tol && it < max.iter) {
            wx <- W %*% X
            gwx <- tanh(alpha * wx)
            v1 <- gwx %*% t(X)/p
            g.wx <- alpha * (1 - (gwx)^2)
            v2 <- Diag(apply(g.wx, 1, FUN = mean)) %*% W
            W1 <- v1 - v2
            sW1 <- La.svd(W1)
            W1 <- sW1$u %*% Diag(1/sW1$d) %*% t(sW1$u) %*% W1
            lim[it + 1] <- max(Mod(Mod(diag(W1 %*% t(W))) - 1))
            W <- W1
            if (verbose)
                message("Iteration ", it, " tol = ", format(lim[it + 1]))
            it <- it + 1
        }
    }
    if (fun == "exp") {
        while (lim[it] > tol && it < max.iter) {
            wx <- W %*% X
            gwx <- wx * exp(-(wx^2)/2)
            v1 <- gwx %*% t(X)/p
            g.wx <- (1 - wx^2) * exp(-(wx^2)/2)
            v2 <- Diag(apply(g.wx, 1, FUN = mean)) %*% W
            W1 <- v1 - v2
            sW1 <- La.svd(W1)
            W1 <- sW1$u %*% Diag(1/sW1$d) %*% t(sW1$u) %*% W1
            lim[it + 1] <- max(Mod(Mod(diag(W1 %*% t(W))) - 1))
            W <- W1
            if (verbose)
                message("Iteration ", it, " tol = ", format(lim[it + 1]))
            it <- it + 1
        }
    }
    return(W)
}
