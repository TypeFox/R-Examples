fastICA <-
function (X, n.comp, alg.typ = c("parallel","deflation"),
          fun = c("logcosh", "exp"),
          alpha = 1, method = c("R", "C"),
          row.norm = FALSE, maxit = 200, tol = 1e-04,
          verbose = FALSE, w.init=NULL)
{
    dd <- dim(X)
    d <- dd[dd != 1L]
    if (length(d) != 2L)
        stop("data must be matrix-conformal")
    X <- if (length(d) != length(dd)) matrix(X, d[1L], d[2L])
    else as.matrix(X)

    if (alpha < 1 || alpha > 2)
        stop("alpha must be in range [1,2]")
    method <- match.arg(method)
    alg.typ <- match.arg(alg.typ)
    fun <- match.arg(fun)
    n <- nrow(X)
    p <- ncol(X)

    if (n.comp > min(n, p)) {
        message("'n.comp' is too large: reset to ", min(n, p))
        n.comp <- min(n, p)
    }
    if(is.null(w.init))
        w.init <- matrix(rnorm(n.comp^2),n.comp,n.comp)
    else {
        if(!is.matrix(w.init) || length(w.init) != (n.comp^2))
            stop("w.init is not a matrix or is the wrong size")
    }
    if (method == "R") {
        if (verbose) message("Centering")

        X <- scale(X, scale = FALSE)

        X <- if (row.norm) t(scale(X, scale=row.norm)) else t(X)

        if (verbose) message("Whitening")
        V <- X %*% t(X)/n

        s <- La.svd(V)
        D <- diag(c(1/sqrt(s$d)))

        K <- D %*% t(s$u)
        K <- matrix(K[1:n.comp, ], n.comp, p)
        X1 <- K %*% X
        a <- if (alg.typ == "deflation")
            ica.R.def(X1, n.comp, tol = tol, fun = fun,
                      alpha = alpha, maxit = maxit, verbose = verbose, w.init = w.init)
        else if (alg.typ == "parallel")
            ica.R.par(X1, n.comp, tol = tol, fun = fun,
                      alpha = alpha, maxit = maxit, verbose = verbose, w.init = w.init)
        w <- a %*% K
        S <- w %*% X
        A <- t(w) %*% solve(w %*% t(w))
        return(list(X = t(X), K = t(K), W = t(a), A = t(A), S = t(S)))
    } else if (method == "C") {
        a <- .C(icainc_JM,
                as.double(X),
                as.double(w.init),
                as.integer(p),
                as.integer(n),
                as.integer(n.comp),
                as.double(alpha),
                as.integer(1),
                as.integer(row.norm),
                as.integer(1L + (fun == "exp")),
                as.integer(maxit),
                as.double(tol),
                as.integer(alg.typ != "parallel"),
                as.integer(verbose),
                X = double(p * n),
                K = double(n.comp * p),
                W = double(n.comp * n.comp),
                A = double(p * n.comp),
                S = double(n.comp * n))
        X1 <- matrix(a$X, n, p)
        K <- matrix(a$K, p, n.comp)
        W <- matrix(a$W, n.comp, n.comp)
        A <- matrix(a$A, n.comp, p)
        S <- matrix(a$S, n, n.comp)
        list(X = X1, K = K, W = W, A = A, S = S)
    }
}

ica.R.def <-
    function (X, n.comp, tol, fun, alpha, maxit, verbose, w.init)
{
    if (verbose && fun == "logcosh")
        message("Deflation FastICA using logcosh approx. to neg-entropy function")
    if (verbose && fun =="exp")
        message("Deflation FastICA using exponential approx. to neg-entropy function")
    n <- nrow(X)
    p <- ncol(X)
    W <- matrix(0, n.comp, n.comp)
    for (i in 1:n.comp) {
        if (verbose) message("Component ", i)
        w <- matrix(w.init[i,], n.comp, 1)
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
        lim <- rep(1000, maxit)
        it <- 1
        if (fun == "logcosh") {
            while (lim[it] > tol && it < maxit) {
                wx <- t(w) %*% X
                gwx <- tanh(alpha * wx)
                gwx <- matrix(gwx, n.comp, p, byrow = TRUE)
                xgwx <- X * gwx
                v1 <- apply(xgwx, 1, FUN = mean)
                g.wx <- alpha * (1 - (tanh(alpha * wx))^2)
                v2 <- mean(g.wx) * w
                w1 <- v1 - v2
                w1 <- matrix(w1, n.comp, 1)
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
                w <- matrix(w1, n.comp, 1)
            }
        }
        if (fun == "exp") {
            while (lim[it] > tol && it < maxit) {
                wx <- t(w) %*% X
                gwx <- wx * exp(-(wx^2)/2)
                gwx <- matrix(gwx, n.comp, p, byrow = TRUE)
                xgwx <- X * gwx
                v1 <- apply(xgwx, 1, FUN = mean)
                g.wx <- (1 - wx^2) * exp(-(wx^2)/2)
                v2 <- mean(g.wx) * w
                w1 <- v1 - v2
                w1 <- matrix(w1, n.comp, 1)
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
                w <- matrix(w1, n.comp, 1)
            }
        }
        W[i, ] <- w
    }
    W
}

ica.R.par <- function (X, n.comp, tol, fun, alpha, maxit, verbose, w.init)
{
    Diag <- function(d) if(length(d) > 1L) diag(d) else as.matrix(d)
    n <- nrow(X)
    p <- ncol(X)
    W <- w.init
    sW <- La.svd(W)
    W <- sW$u %*% Diag(1/sW$d) %*% t(sW$u) %*% W
    W1 <- W
    lim <- rep(1000, maxit)
    it <- 1
    if (fun == "logcosh") {
        if (verbose)
            message("Symmetric FastICA using logcosh approx. to neg-entropy function")
        while (lim[it] > tol && it < maxit) {
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
        if (verbose)
            message("Symmetric FastICA using exponential approx. to neg-entropy function")
        while (lim[it] > tol && it < maxit) {
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
    W
}
