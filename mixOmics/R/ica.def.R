# The function ica.par and ica.def are borrowed from the fastICA package (see references).

ica.def <-
    function (X, ncomp, tol, fun, alpha, max.iter, verbose, w.init)
{
    n <- nrow(X)
    p <- ncol(X)
    W <- matrix(0, ncomp, ncomp)
    for (i in 1:ncomp) {
        w <- matrix(w.init[i,], ncomp, 1)
        if (i > 1) {
            t <- w
            t[1:length(t)] <- 0
            for (u in 1:(i - 1)) {
                k <- sum(w * W[u, ])
                t <- t + k * W[u, ]
            }
            w <- w - t
        }
        w <- w/sqrt(sum(w^2))
        lim <- rep(1000, max.iter)
        it <- 1
        if (fun == "logcosh") {
            while (lim[it] > tol && it < max.iter) {
                wx <- t(w) %*% X
                gwx <- tanh(alpha * wx)
                gwx <- matrix(gwx, ncomp, p, byrow = TRUE)
                xgwx <- X * gwx
                v1 <- apply(xgwx, 1, FUN = mean)
                g.wx <- alpha * (1 - (tanh(alpha * wx))^2)
                v2 <- mean(g.wx) * w
                w1 <- v1 - v2
                w1 <- matrix(w1, ncomp, 1)
                it <- it + 1
                if (i > 1) {
                    t <- w1
                    t[1:length(t)] <- 0
                    for (u in 1:(i - 1)) {
                        k <- sum(w1 * W[u, ])
                        t <- t + k * W[u, ]
                    }
                    w1 <- w1 - t
                }
                w1 <- w1/sqrt(sum(w1^2))
                lim[it] <- Mod(Mod(sum((w1 * w))) - 1)
                w <- matrix(w1, ncomp, 1)
            }
        }
        if (fun == "exp") {
            while (lim[it] > tol && it < max.iter) {
                wx <- t(w) %*% X
                gwx <- wx * exp(-(wx^2)/2)
                gwx <- matrix(gwx, ncomp, p, byrow = TRUE)
                xgwx <- X * gwx
                v1 <- apply(xgwx, 1, FUN = mean)
                g.wx <- (1 - wx^2) * exp(-(wx^2)/2)
                v2 <- mean(g.wx) * w
                w1 <- v1 - v2
                w1 <- matrix(w1, ncomp, 1)
                it <- it + 1
                if (i > 1) {
                    t <- w1
                    t[1:length(t)] <- 0
                    for (u in 1:(i - 1)) {
                        k <- sum(w1 * W[u, ])
                        t <- t + k * W[u, ]
                    }
                    w1 <- w1 - t
                }
                w1 <- w1/sqrt(sum(w1^2))
                lim[it] <- Mod(Mod(sum((w1 * w))) - 1)
                if (verbose)
                    message("Iteration ", it - 1, " tol = ", format(lim[it]))
                w <- matrix(w1, ncomp, 1)
            }
        }
        W[i, ] <- w
    }
    return(W)
}
